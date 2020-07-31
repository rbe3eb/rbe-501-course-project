function exampleCommandMoveToTaskConfig(coordinator, taskConfig, tspan)
        syms tt q1 q2 q3 q4 q5 q6 q7 real;
        
        ti = tspan(1);
        tf = tspan(2);
        
        HT = coordinator.HT;
        Teg = coordinator.Teg;
        Jv = coordinator.Jv;
        Jw = coordinator.Jw;
        MotionModel = coordinator.MotionModel;

        % Current robot joint configuration
        jointInit = coordinator.CurrentRobotJConfig;
        currentRobotJConfig = wrapToPi(jointInit');
        
        % Final (desired) end-effector pose
        currP = eval(subs(HT{end}(1:3,4),[q1,q2,q3,q4,q5,q6,q7],currentRobotJConfig'));
        currR = eval(subs(HT{end}(1:3,1:3),[q1,q2,q3,q4,q5,q6,q7],currentRobotJConfig'));
        
        %taskConfig = taskConfig*inv(Teg);
        desR = (taskConfig(1:3,1:3));
        desP = taskConfig(1:3,4);
        
        jointFinal = IK(HT{end},Jv{end},Jw{end},desP,desR,currP,currR,jointInit);
        
        qi = jointInit'; qf = jointFinal;
        dqi = zeros(7,1); dqf = zeros(7,1);
        
        [qEqn, dqEqn, ddqEqn] = jointSpaceTrajectory(HT{end},Jv{end},qi,dqi,qf,dqf,ti,tf);
        
        initState = [qi; dqi];
        
        % Simulate motion jointSpaceMotionControl(qEqn,dqEqn,ddqEqn,HT,Jv,qi,dqi,ti,tf) 
        options = odeset('RelTol', 1e-2, 'AbsTol', [1e-2*ones(14,1)]);
        [T,X] = ode15s(@(t,x) exampleHelperTimeBasedStateInputsPickPlace(t,x,MotionModel,qEqn,dqEqn,ddqEqn,HT,Teg,Jv), [ti tf], initState,options);
        
        % Forward kinematics for end-effector and gripper
        cnt = size(T,1);
        eeT = cell(cnt,1); gripperT = cell(cnt,1);
        Xdes = zeros(cnt,14);
        eePos = zeros(cnt,3); eePosDes = zeros(cnt,3);
        for n=1:cnt
            t = T(n);
            q = X(n,1:7);
            dq = X(n,8:14);
            
            Xdes(n,:) = [eval(subs(qEqn,tt,t))', eval(subs(dqEqn,tt,t))'];
            
            eePos(n,:) = eval(subs(HT{end}(1:3,4),[q1,q2,q3,q4,q5,q6,q7],q));
            eePosDes(n,:) = eval(subs(HT{end}(1:3,4),[q1,q2,q3,q4,q5,q6,q7],Xdes(n,1:7)));     
        end
        
        plotData(T,X,Xdes,eePos,eePosDes,'PD Controller','Joint Space');

        % Uncomment below to display all successful outputs
        % disp('Executing collision-free trajectory...')
        visualizePath(coordinator,eeT);
        visualizeRobot(coordinator, X(:,1:7), eeT, T);
        
        % Deleta path on plot
        coordinator.PathHandle.Visible = 'off';

        % Update current robot configuration
        coordinator.CurrentRobotJConfig = X(end,1:7);
        coordinator.CurrentRobotTaskConfig = eeT{end}; 

        % Trigger Stateflow chart Event
        %coordinator.FlowChart.taskConfigReached;
end