function exampleCommandMoveToTaskConfig(coordinator, taskConfig, tolerance, avoidCollisions)
        syms q1 q2 q3 q4 q5 q6 q7 real;
        
        HT = coordinator.HT;
        Teg = coordinator.Teg;
        Jv = coordinator.Jv;
        Jw = coordinator.Jw;
        MotionModel = coordinator.MotionModel;
        

        % Current robot joint configuration
        jointInit = coordinator.CurrentRobotJConfig;
        currentRobotJConfig = wrapToPi(jointInit');
        
        % Final (desired) end-effector pose
        currP = subs(HT{end}(1:3,4),[q1,q2,q3,q4,q5,q6,q7],currentRobotJConfig');
        currR = subs(HT{end}(1:3,1:3),[q1,q2,q3,q4,q5,q6,q7],currentRobotJConfig');
        
        desR = (taskConfig(1:3,1:3));
        desP = taskConfig(1:3,4);
        jointFinal = IK(HT{end},Jv{end},Jw{end},desP,desR,currP,currR,jointInit);
        
        qi = jointInit'; qf = jointFinal;
        dqi = zeros(7,1); dqf = zeros(7,1);
        ti = 0; tf = 5;
        
        [qEqn, dqEqn, ddqEqn] = jointSpaceTrajectory(HT{end},Jv{end},qi,dqi,qf,dqf,ti,tf);
        
        initState = [qi; dqi];
        
        % Simulate motion jointSpaceMotionControl(qEqn,dqEqn,ddqEqn,HT,Jv,qi,dqi,ti,tf) 
        [T,X] = ode15s(@(t,x) exampleHelperTimeBasedStateInputsPickPlace(t,x,MotionModel,qEqn,dqEqn,ddqEqn,HT,Teg,Jv), [ti tf], initState);
        
        % Forward kinematics for end-effector and gripper
        cnt = size(T,1);
        eeT = cell(cnt,1); gripperT = cell(cnt,1);
        for n=1:cnt
            q = X(n,1:7);
            eeT{n} = eval(subs(HT{end},[q1,q2,q3,q4,q5,q6,q7],q));
            gripperT{n} = eeT{n}*Teg;
        end

        % Uncomment below to display all successful outputs
        % disp('Executing collision-free trajectory...')
        visualizePath(coordinator,gripperT);
        visualizeRobot(coordinator, X(:,1:7), eeT, T);
        
        % Deleta path on plot
        coordinator.PathHandle.Visible = 'off';

        % Update current robot configuration
        coordinator.CurrentRobotJConfig = X(end,1:7);
        coordinator.CurrentRobotTaskConfig = eeT{end}; 

        % Trigger Stateflow chart Event
        coordinator.FlowChart.taskConfigReached; 
end