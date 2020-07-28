classdef exampleHelperManipCollsionsPickPlace < robotics.manip.internal.InternalAccess
    % Copyright 2020 The MathWorks, Inc.
    % This file is for internal use only and may be modified or removed in
    % a future release.
    %
    %exampleHelperManipCollsions Tools for manipulator collisions
    %   This example helper has tools for collision checking. The object is
    %   constructed as a collision tree using either primitives or meshes.
    %   Helper methods are used for self-collisions and environmental
    %   collisions.
    %
   
    
    %   Copyright 2020 The MathWorks, Inc.
    
    properties
        SourceRigidBodyTree
        
        VizRigidBodyTree
        
        RigidBodyCollisionArray
        
        PartMask
        
        %ExhaustiveChecking   - Boolean that indicates when to stop checking
        %   When the parameter is set to FALSE, the collision-checking
        %   methods exit as soon as a collision is found. When this
        %   parameter is set to TRUE, collision-checking continues between
        %   all applicable bodies, rather than exiting at the first
        %   collision.
        ExhaustiveChecking = false
    end
    
    %% Constructor methods
    
    methods
        function obj = exampleHelperManipCollsionsPickPlace(tree, partMask)
            
            % Define the tree from the input
            validateattributes(tree, {'robotics.RigidBodyTree'}, {'scalar'}, 'exampleHelperManipCollsions', 'tree');
            obj.SourceRigidBodyTree = copy(tree);
            obj.VizRigidBodyTree = copy(tree);

            % Initialize collision array from tree
            obj.RigidBodyCollisionArray = cell(tree.NumBodies+1, 2);

            % Update the data format
            obj.SourceRigidBodyTree.DataFormat = 'column';
            obj.VizRigidBodyTree.DataFormat = 'column';
                        
            % Create the collision object array using the visuals
            % in the associated rigid body tree
            obj.createCollisionArrayFromVisuals; 
            
            obj.PartMask  = partMask;
        end
    end
    
    %% Public collision-checking methods
    
    methods        
        % Skip self-collision
        
         function skipCollisionCheck(obj, bodyIndex)
            %skipCollisionCheck
            obj.PartMask(bodyIndex) = 0;
         end
         
         function enableCollisionCheck(obj, bodyIndex)
            %enableCollisionCheck
            obj.PartMask(bodyIndex) = 1;
         end
        
        %checkRobotWorldCollision Check collision with world
        function [allDistances, allBodiesWtnPts, allBodyNames] = checkRobotWorldCollision(obj, config, worldCollisionObjects)
            
             
            % Populate transformTree, which is a cell array of all body
            % transforms w.r.t. global frame
            transformTree = obj.getTransformTreeInternal(config);
            
            % Pairwise checking
            numBodies = size(obj.RigidBodyCollisionArray,1);
            nonBaseBodies = numBodies-1;
            actualBodiesNumber = numel(nonzeros(obj.PartMask));
            allBodyNames = cell(1, actualBodiesNumber);
            
            numObstacles = numel(worldCollisionObjects);
            allDistances = zeros(actualBodiesNumber*numObstacles,1);
            allBodiesWtnPts = zeros(3, 2, actualBodiesNumber, numObstacles );
            distIter = 1;
            bodiesiter = 1;
            
            for i = 1:nonBaseBodies                
                % Iterate through all robot bodies. If no geometry is
                % assigned, skip this body
                if ~isempty(obj.RigidBodyCollisionArray{i+1}) && obj.PartMask(i)==1
                    obj.RigidBodyCollisionArray{i+1,1}.Pose = transformTree{i+1};
                    
                    % Iterate through all collision geometries
                    for j = 1:numel(worldCollisionObjects)
                        [localCollisionStatus, sepDist, wPts] = checkCollision(obj.RigidBodyCollisionArray{i+1,1}, worldCollisionObjects{j});
                        if localCollisionStatus
                            allDistances(distIter) = 0;
                        else
                            allBodiesWtnPts(:,:,i,j) = wPts;
                            allDistances(distIter) = sepDist;
                        end
                        distIter = distIter + 1;
                    end
                allBodyNames{bodiesiter} = obj.SourceRigidBodyTree.BodyNames{i}; 
                bodiesiter = bodiesiter + 1; 
                end
            end
        end
    end
    %% Helper Methods
    
    methods (Access = private)
        function createCollisionArrayFromVisuals(obj)
            %createCollisionArrayFromVisuals
            
            % Use the source rigid body tree to get basic info
            tree = obj.SourceRigidBodyTree;
            
            simplificationThreshold = 50;
            % For each of the bodies, get the body internal, which links to
            % the visual information
            allBodies = [{tree.Base} tree.Bodies];
            for i = 1:numel(allBodies)
                if i == 1
                    bodyInternal = tree.Base.BodyInternal;
                else
                    bodyInternal = tree.Bodies{i-1}.BodyInternal;
                end
                
                visualsInternal = bodyInternal.VisualsInternal;
                if ~isempty(visualsInternal)
                
                    V = visualsInternal{1}.Vertices; % only extracting the first geometry in each body
                    T = visualsInternal{1}.Tform;
                    V = double(V);

                    originalNumV = size(V,1);
                    VT = T*[V'; ones(1,originalNumV)];
                    V = VT(1:3,:)';
                    V2 = V;

                    if originalNumV > simplificationThreshold
                        %simplify the mesh
                        f = figure;
                        f.Visible = "off";
                        ax = axes(f);
                        F = convhull(V(:,1), V(:,2), V(:,3));
                        V1 = V(unique(F(:)),:);
                        F1 = convhull(V1(:,1), V1(:,2), V1(:,3));
                        if size(V1,1) <= simplificationThreshold
                            V2 = V1;
                        else
                            p = patch(ax, 'Faces',F1,'Vertices',V1);
                            [~, V2] = reducepatch(p, simplificationThreshold/size(V1,1));
                        end
                        close(f)
                    else
                        V2 = V;
                    end

                    obj.RigidBodyCollisionArray{i,1} = collisionMesh(V2);
                    obj.RigidBodyCollisionArray{i,2} = T;
                end
            end
        end
        
     
        function assignCollisionPose(obj, config)
            %assignCollisionPose Assign poses to collision objects
            %   This method updates the cell array of collision objects
            %   with the pose of the rigid body tree given the input
            %   specified by the CONFIG vector. This method uses the
            %   internal transform tree method, which is essentially an
            %   efficient variation on iterative calls to getTransform, and
            %   then uses that to update the Pose property of each
            %   corresponding collision object.
            
            transformTree = obj.getTransformTreeInternal(config);
            for i = 1:size(obj.RigidBodyCollisionArray,1)
                if ~isempty(obj.RigidBodyCollisionArray{i,1})
                    % The pose is the product of the body frame's transform
                    % relative to the base, and the collision object's
                    % transform relative to the body frame.
                    jointToCollisionMeshTransform = obj.RigidBodyCollisionArray{i,2};
                    obj.RigidBodyCollisionArray{i,1}.Pose = transformTree{i}*jointToCollisionMeshTransform;
                end
            end
        end
        
        function transformTree = getTransformTree(obj, config)
            %getTransformTree Compute position of each body in tree
            %   This method outputs a cell array of the transforms of each
            %   body in the associated rigidBodyTree. This method is
            %   inefficient, and an internal method below avoids redundant
            %   calls.
            
            % Initialize for the base
            transformTree = {eye(4)};
            
            % Fill in the other body transforms
            for i = 1:obj.SourceRigidBodyTree.NumBodies
                tree = obj.SourceRigidBodyTree;
                transformTree{i+1} = ... 
                    getTransform(tree, config, tree.Base.Name, tree.Bodies{i}.Name); % tree.Bodies{i} is very inefficient
            end
        end
        
        function Ttree = getTransformTreeInternal(obj, config)
            %getTransformTree Compute position of each body in tree
            %   This is an internal helper method that efficiently computes
            %   the positions of all the bodies in the rigid body tree.
            %   This is faster than computing getTransform in a loop, as
            %   that method calls each of the transforms for each body
            %   leading up to that chain for every call (so it ends up
            %   with a bunch of redundant calls).
            
            % Use an internal method and concatentate with the known
            % transform of the base, which is at the origin
            Ttree = [{eye(4)} obj.SourceRigidBodyTree.TreeInternal.forwardKinematics(config)];
        end
    end
end

