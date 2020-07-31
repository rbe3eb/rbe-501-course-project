clear all; close all; clc;

load('exampleHelperKINOVAGen3RobotiqGripper.mat'); 

coordinator = init(robot);

% Pick and place parts
while true
    part = pickObject(coordinator);
    placeObject(coordinator,part);
end


function [coordinator] = init(robot)
    currentRobotJConfig = homeConfiguration(robot);
    coordinator = exampleHelperCoordinatorPickPlace(robot,currentRobotJConfig, "gripper");
    %go to home configuration
    coordinator.HomeRobotTaskConfig = trvec2tform([0.4, 0, 0.6])*axang2tform([0 1 0 pi]);
    %hold over white beaker to pour
    coordinator.HoldConfig = trvec2tform([[0.4 0.15 0.50]])*axang2tform([0 1 0 pi]);
    %place blue beaker back on table
    coordinator.PlacingPose{1} = trvec2tform([[-0.05, 0.45, 0.32]])*axang2tform([0 1 0 pi]);
    %place red beaker back on table
    coordinator.PlacingPose{2} = trvec2tform([[-0.25, 0.45, 0.32]])*axang2tform([0 1 0 pi]);
    %place white beaker on mixing platform shelf
    coordinator.PlacingPose{3} = trvec2tform([[0.3 0.5 0.345]])*axang2tform([0 1 0 pi]);
    %place white beaker back in center table area
    %coordinator.PlacingPose{4} = trvec2tform([[0.4 0.15 0.32]])*axang2tform([0 1 0 pi]);
    % Build world
    exampleCommandBuildWorld(coordinator);
    % Move to home position
    exampleCommandMoveToTaskConfig(coordinator, coordinator.HomeRobotTaskConfig, [0 3]);
    % Detect parts
    exampleCommandDetectParts(coordinator);
    % Classify parts
    exampleCommandClassifyParts(coordinator);
end

function [] = placeObject(coordinator,part)
    fprintf('\nPlacing Object\n')
    % Move to placement approach position
% Move to placement approach position
    %exampleCommandMoveToTaskConfig(coordinator, coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}*trvec2tform([0, 0, -0.2]), [0 1]);

    % Move to placement approach position
    %exampleCommandMoveToTaskConfig(coordinator, coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}, [0 1]);
    currT = coordinator.CurrentRobotTaskConfig;
    desT = coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt-[0;0;0.5]};
    moveObject(coordinator, currT, desT);

    % Deactivate gripper
    exampleCommandActivateGripper(coordinator,'off')

    % Move to retracted position
    exampleCommandMoveToTaskConfig(coordinator, coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}*trvec2tform([0, 0, -0.2]), [0 1]);

    % Update the placing pose so that the next part is placed elsewhere
    coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}(1,4) =  coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}(1,4) - 0.15;
    fprintf('Part has been placed\n')
end

function [currPart] = pickObject(coordinator)
        syms tt
        currPart = coordinator.NextPart+1;
        if currPart > 3
            currPart = 3;
        end
        updateFlow(currPart);
        exampleCommandPickingLogic(coordinator)
        % Compute grasp pose
        exampleCommandComputeGraspPose(coordinator);
        % Move to picking approach position
        
        %exampleCommandMoveToTaskConfig(coordinator, coordinator.GraspPose*trvec2tform([0,0,-0.1]),[0 2]);
        % Move to reach position
        currT = coordinator.CurrentRobotTaskConfig;
        graspT = coordinator.GraspPose;
        moveObject(coordinator, currT, graspT);
        fprintf('\nObject reached\n')

        % Activate the gripper
        exampleCommandActivateGripper(coordinator,'on');
        % Move to retracted position
        exampleCommandMoveToTaskConfig(coordinator,coordinator.GraspPose*trvec2tform([0, 0, -0.2]), [0 1]);
        fprintf('\nObject picked up\n')
        % Move to retracted position
        currT = coordinator.CurrentRobotTaskConfig;
        graspT = coordinator.HoldConfig;
        moveObject(coordinator, currT, graspT);
        updateFlow(currPart*-1);
end

function updateFlow(part)
    if part==1
        fprintf("\nPicking up distilled water\n")
    elseif part==2
        fprintf("\nPicking up NaCL\n")
    elseif part==3
        fprintf("\nPicking up distilled Main Beaker\n")
    elseif part==-1
        fprintf("\nPouring water into beaker...\n")
    elseif part==-2
        fprintf("\nPouring NaCL into Beaker\n")
    end
end

function moveObject(coordinator, currT, graspT)
        syms tt
        toolSpeed = 0.25; % m/s
        distance = norm(tform2trvec(currT)-tform2trvec(graspT));
        tf = ceil(distance/toolSpeed); ti = 0;
        Xi = currT(1:3,4); Xf = graspT(1:3,4);
        dXi = zeros(3,1); dXf = zeros(3,1);
        [pEqn,~,~] = taskSpaceTrajectory(Xi,Xf,dXi,dXf,ti,tf);
        for t=linspace(ti+0.2,tf,tf*1.5)
            desP = eval(subs(pEqn,tt,t));
            refPose = [graspT(1:4,1:3), [desP; 1]];
            exampleCommandMoveToTaskConfig(coordinator, refPose, [t t+1]);     
        end
end