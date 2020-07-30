close all; clear all; clc;

load('exampleHelperKINOVAGen3RobotiqGripper.mat'); 

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

coordinator.FlowChart = exampleHelperFlowChartPickPlace('coordinator', coordinator); 

% Trigger event to start Pick and Place in the Stateflow Chart
coordinator.FlowChart.startPickPlace; 