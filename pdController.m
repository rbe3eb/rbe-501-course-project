function [] = pdController()
    clc; clear; close all;
    global theta theta_d desPos desTheta desTheta_d T X;
    theta = []; theta_d = []; desPos = []; desTheta = []; desTheta_d = [];
    torque = []; eePos = []; eePosErr = []; joint2Pos = []; tf = 10;
    
    fprintf('Loading Variables\n')
    [T,M,C,G,q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6] = loadVars();
    fprintf('Variables loaded!\n')
    
    setPointTrajectory();

    function dx = planarArmODE(t,x)
        syms q1 q2 q3 q4 q5 q6 ...
            dq1 dq2 dq3 dq4 dq5 dq6 real
        q = [q1,q2,q3,q4,q5,q6];
        dq = [dq1,dq2,dq3,dq4,dq5,dq6];
        
        theta = x(1:6,1); % Current Position
        theta_d = x(7:12,1); % Current Velocity
        M = eval(subs(M, [q,dq], [theta' theta_d']));
        C = eval(subs(C, [q,dq], [theta' theta_d']));
        G = eval(subs(G, [q,dq], [theta' theta_d']));
        % invert M here to save time
        invM = inv(M);
        % PD controller
        tau = PDControl(theta, theta_d);
        % Append Torque
        torque = [torque, tau];
        dx = zeros(6,1);
        dx(1:6) = x(7:12);
        dx(7:12) = (invM*tau) - ((invM*C).*x(7:12));
        % Compute end-effector position
        eePos = [eePos, computeFK(T{end},theta)]
        joint2Pos = [joint2Pos, computeFK(T{end-1},theta)];
        eePosErr = [eePosErr, desPos-eePos(1:3,end)];
    end
    
    % Controller
    function tau = PDControl(theta, theta_d)
        Kp = 10*ones(6,6);
        Kv = 10*ones(6,6);
        e = desTheta - theta; % position error
        de = desTheta_d - theta_d; % velocity error
        tau = Kp*e + Kv*de;
    end
    
    % Given cartesian coordinates, compute joint angles
    function q = computeIK(pos)
        syms q1 q2 q3 q4 q5 q6 real
        q = NaN*ones(6,1);
        eqn = T{end}(1:3,4) == [pos(1); pos(2); pos(3)];
        S = solve(eqn, [q1,q2,q3,q4,q5,q6]);
        if size(S.q1,1) == 0
            return;
        end
        q(1) = S.q1(1); q(2) = S.q2(1); q(3) = S.q3(1);
        q(4) = S.q4(1); q(5) = S.q5(1); q(6) = S.q6(1);
    end
    
    % Given joint angles, return cartesian coordinates
    function pos = computeFK(TT,q)
        pos = (subs(TT(1:3,4), [q1,q2,q3,q4,q5,q6], [q(1),q(2),q(3) ...
            q(4),q(5),q(6)]));
    end

    function plotData()
        figure('Name','Joint Angles vs Time under PD Control, GComp');
        for n=1:6
            plot(T, X(:,n),'r-');
            hold on
        end
        xlim([0, tf]);
        xlabel('time (secs)');
        ylabel('joint angles (rads)');
        legend('q1', 'q2', 'q3', 'q4', 'q5', 'q6')
        title('Joint Angles vs Time under PD Control, GComp');
        hold off
        drawnow
        
        figure('Name','Joint Velocities vs Time under PD Control,GComp');
        for n=7:12
           plot(T, X(:,n),'r-');
            hold on 
        end
        xlim([0, tf]);
        xlabel('time (secs)')
        ylabel('joint velocities (rads/secs)')
        legend('q1_dot', 'q2_dot', 'q3_dot','q4_dot', 'q5_dot', 'q6_dot')
        title('Joint Velocities vs Time under PD Control, GComp');
        hold off
        drawnow
        
        figure('Name','Control Input vs Time under PD control, GComp');
        for n=1:6
        	plot(T, torque(n,1:size(T,1)),'r-' );
            hold on
        end
        xlim([0, tf]);
        xlabel('time (secs)')
        ylabel('control input (N-m)')
        legend('torque1', 'torque2', 'torque3','torque4', 'torque5', ...
            'torque6')
        title('Control Input vs Time under PD control, GComp');
        hold off
        drawnow

        figure('Name','EE Position under PD control, GComp');
        for n=1:6
            plot(T, eePos(n,1:size(T,1)), 'r-');
            hold on 
        end
        xlim([0, tf]);
        xlabel('time (secs)')
        ylabel('end-effector positions (m)')
        legend('X position', 'Y position', 'Z position')
        title('EE Position vs Time under PD control, GComp')
        hold off
        drawnow
        
        figure('Name','EE Position Error vs Time under PD control, GComp');
        for n=1:6
            plot(T, eePosErr(n,1:size(T,1)), 'r-');
            hold on
        end
        xlim([0, tf]);
        xlabel('time (secs)')
        ylabel('error in ee positions (m)')
        legend('X position error', 'Y position error', 'Z position error')
        title('EE Position Error vs Time under PD control, GComp')
        hold off
        drawnow
        
        %Set up stick model
        end_eff_x = eePos(1,:);
        end_eff_y = eePos(2,:);
        end_eff_z = eePos(3,:);
        joint1_x = 0;
        joint1_y = 0;
        joint1_z = 0.3;
        joint2_x = joint2Pos(1,:);
        joint2_y = joint2Pos(2,:);
        joint2_z = joint2Pos(3,:);
        line1_x = [0 joint1_x(1,1) joint2_x(1,1) end_eff_x(1,1)];
        line1_y = [0 joint1_y(1,1) joint2_y(1,1) end_eff_y(1,1)];
        line1_z = [0 joint1_z(1,1) joint2_z(1,1) end_eff_z(1,1)];
        line2_x = [0 joint1_x(end,1) joint2_x(end,1) end_eff_x(end,1)];
        line2_y = [0 joint1_y(end,1) joint2_y(end,1) end_eff_y(end,1)];
        line2_z = [0 joint1_z(end,1) joint2_z(end,1) end_eff_z(end,1)];
        
        figure('Name', 'End Effector Trajectory w/ Stick Model')
        plot3(end_eff_x, end_eff_y, end_eff_z, joint2_x, joint2_y, joint2_z, 'r--', line1_x, ...
        line1_y, line1_z, line2_x, line2_y, line2_z);
        ylim([0,tf]);
        xlabel('X (meters)');
        ylabel('Y (meters)');
        zlabel('Z (meters)');
        legend('End effector trajectory', 'Joint 2 trajectory', 'Initial Position' ...
            ,'Final position')
        title('End Effector Trajectory w/ Stick Model Gcomp');
    end
    
    % Compute trajectory from [0, 0.3, 0.6] to [0.3, 0, 0.6]
    function setPointTrajectory()
        % Desired state variables
        desPos = [0; 0.3; 0.6];
        desTheta = computeIK(desPos'); % Desired Set-Point Position
        desTheta_d = [0; 0; 0; 0; 0; 0]; % Desired velocity (Derivative of theta_d)

        % initial state variables
        initPos = [0.3; 0; 0.6];
        initTheta = computeIK(initPos');
        initTheta_d = zeros(3,1);
        initState = [initTheta ; initTheta_d]';
        if (any(isnan([initTheta; desTheta])))
            fprintf('\nThe requested cartesian point: [%0.3f, %0.3f, %0.3f] is not reachable\nExiting...\n',desX,desY,desZ);
            return;
        end
        fprintf('\nComputing trajectory from [%0.3f, %0.3f, %0.3f] to [%0.3f, %0.3f, %0.3f] \n',0,0.3,0.6,0.3,0,0.6);
        [T,X] = ode45(@(t,x)planarArmODE(t,x),[0 tf],initState);
        plotData();
    end
end

function [T,M,C,G,q1,q2,q3,q4,q5,q6,dq1,dq2,dq3,dq4,dq5,dq6] = loadVars()
    syms q1 q2 q3 q4 q5 q5 q6 ...
        dq1 dq2 dq3 dq4 dq5 dq6 real;
    filenames = {'T.mat','M.mat','C.mat','G.mat'};  
    for var = 1:numel(filenames)
        load(filenames{var})
    end 
end