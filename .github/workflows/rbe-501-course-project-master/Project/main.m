close all; clc;

load('exampleHelperKINOVAGen3RobotiqGripper.mat'); 

coordinator = init(robot);

% Pick and place parts
while true
    pickObject(coordinator);
    placeObject(coordinator);
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
    exampleCommandMoveToTaskConfig(coordinator, coordinator.HomeRobotTaskConfig, [0 2]);
    % Detect parts
    exampleCommandDetectParts(coordinator);
    % Classify parts
    exampleCommandClassifyParts(coordinator);
end

function [] = placeObject(coordinator)
    % Move to placement approach position
    exampleCommandMoveToTaskConfig(coordinator, coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}*trvec2tform([0, 0, -0.2]),0.1, true);

    % Move to placement approach position
    exampleCommandMoveToTaskConfig(coordinator, coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt},0.01, false);

    % Deactivate gripper
    exampleCommandActivateGripper(coordinator,'off')

    % Move to retracted position
    exampleCommandMoveToTaskConfig(coordinator, coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}*trvec2tform([0, 0, -0.2]),0.07, false);

    % Update the placing pose so that the next part is placed elsewhere
    coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}(1,4) =  coordinator.PlacingPose{coordinator.DetectedParts{coordinator.NextPart}.placingBelt}(1,4) - 0.15;

    fprintf('Part has been placed\n')
end

function [] = pickObject(coordinator)
        exampleCommandPickingLogic(coordinator)
        % Compute grasp pose
        exampleCommandComputeGraspPose(coordinator);
        % Move to picking approach position
        
        exampleCommandMoveToTaskConfig(coordinator, coordinator.GraspPose*trvec2tform([0,0,-0.1]),[2 3]);
        % Move to reach position
        exampleCommandMoveToTaskConfig(coordinator, coordinator.GraspPose, [3 4]);

        % Activate the gripper
        exampleCommandActivateGripper(coordinator,'on');

        % Move to retracted position
        exampleCommandMoveToTaskConfig(coordinator,coordinator.GraspPose*trvec2tform([0, 0, -0.2]), [4 5]);

        % Move to retracted position
        exampleCommandMoveToTaskConfig(coordinator, coordinator.HoldConfig, [5 6]);
end