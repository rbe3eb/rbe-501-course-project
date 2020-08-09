function exampleCommandActivateGripper(coordinator, state)% This class is for internal use and may be removed in a future release
%
%CommandActivateGripper Command function to activate gripper  
%   This command activates the gripper. However, this action has just two
%   impacts:
%      - The collision array must be updated to include the part that has
%      been picked up, so that trajectory planning includes the part in the
%      obstacle avoidance stage. This step is handled by the "doit" method
%      in this class
%      - The visualization has to be updated to move the object being
%      picked up by the gripper. This is handled by the dispatcher in the
%      visualization
%

% Copyright 2020 The MathWorks, Inc.

        %   Based on the state, decide whether to activate or
        %   deactivate the gripper
       if strcmp(state,'on') == 1
           % Activate gripper 
           coordinator.PartOnRobot = coordinator.NextPart;
       else
           % Deactivate gripper 
           coordinator.PartOnRobot = 0;
       end
       
       % Trigger Stateflow chart Event
       coordinator.FlowChart.nextAction; 
end