function exampleCommandBuildWorld(coordinator)

%CommandBuildWorld Construct the world used for visualization and collision detection
%   This command constructs the environment, consisting of the workstation,
%   belts, parts to be moved, and obstacles. The robot is constructed
%   separately. These objects are created using collision primitives and
%   placed in the right spots. In the workflow, this world is used during
%   trajectory generation/optimization for collision-checking, and during
%   the visualization step.
%

% Copyright 2020 The MathWorks, Inc.

    % Construct the Workstation (only for visualization)
    bench = collisionBox(0.5, 0.6, 0.05);
    belt1 = collisionBox(1.3, 0.4, 0.05);

    % Place the bench and belts in the proper locations
    TBench = trvec2tform([0.4 0.1 0.2]);
    TBelt1 = trvec2tform([0 0.6 0.2]);

    bench.Pose = TBench;
    belt1.Pose = TBelt1;

    coordinator.World = {bench, belt1};

    % Parts Sizes
    %Main mixing beaker
    box2 = collisionCylinder(0.07, 0.1);
    %Top Shelf
    box1 = collisionCylinder(0.04, 0.07);
    box5 = collisionCylinder(0.04, 0.07);
    %Middle Shelf
    box3 = collisionCylinder(0.04, 0.07);
    box4 = collisionCylinder(0.04, 0.07);
    %Bottom shelf
    box6 = collisionCylinder(0.04, 0.07);
    box7 = collisionCylinder(0.04, 0.07);


    % Parts Positions
    %Main mixing beaker
    TBox2 = trvec2tform([0.4 0.15 0.28]);
    %Top Shelf
    TBox1 = trvec2tform([-0.05 0.55 0.62]);
    TBox5 = trvec2tform([-0.25 0.55 0.62]);
    %Middle Shelf
    TBox3 = trvec2tform([-0.25 0.55 0.44]);
    TBox4 = trvec2tform([-0.05 0.55 0.44]);
    %Bottom shelf
    TBox6 = trvec2tform([-0.05 0.55 0.25]);
    TBox7 = trvec2tform([-0.25 0.55 0.25]);

    box1.Pose = TBox1;
    box2.Pose = TBox2;
    box3.Pose = TBox3;
    box4.Pose = TBox4;
    box5.Pose = TBox5;
    box6.Pose = TBox6;
    box7.Pose = TBox7;

    % Set the part mesh and color
    part1.mesh = box2;
    part2.mesh = box3;
    part3.mesh = box1;
    part4.mesh = box4;
    part5.mesh = box5;
    part6.mesh = box6;
    part7.mesh = box7;

    part1.color = 'white';
    part2.color = 'r';
    part3.color = 'g';
    part4.color = 'blue';
    part5.color = 'orange';
    part6.color = 'purple';
    part7.color = 'y';

    part1.centerPoint = tform2trvec(part1.mesh.Pose);
    part2.centerPoint = tform2trvec(part2.mesh.Pose);
    part3.centerPoint = tform2trvec(part3.mesh.Pose);
    part4.centerPoint = tform2trvec(part4.mesh.Pose);
    part5.centerPoint = tform2trvec(part5.mesh.Pose);
    part6.centerPoint = tform2trvec(part6.mesh.Pose);
    part7.centerPoint = tform2trvec(part7.mesh.Pose);

    part1.plot = [];
    part2.plot = [];
    part3.plot = [];
    part4.plot = [];
    part5.plot = [];
    part6.plot = [];
    part7.plot = [];

    coordinator.Parts = {part1, part2, part3, part4, part5, part6, part7};

    % Visualize world and parts
    visualizeWorld(coordinator)
    visualizeParts(coordinator)

   % Trigger Stateflow chart Event
   coordinator.FlowChart.worldBuilt;
end

