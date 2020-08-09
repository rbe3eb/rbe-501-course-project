function [] = problem_1()
    clear all; close all; clc;
    syms q1 q2 q3 real;
    
    % Forward kinematics
    [HT,linkPos] = hTrans();
    
    % Jacobian
    Jv = Jacobian(linkPos);
    
    % Trajectory beginning and end states
    Xi = [150; 150; 150]; Xf = [150; 150; 0];
    dXi = [0; 0; 0]; dXf = [0; 0; 0];
    tspan = [0 10];
    
    % Generate the task space trajectory
    [pEqn,dpEqn,ddpEqn] = taskSpaceTrajectory(Xi,Xf,dXi,dXf,tspan);
    
    % Start robot at in its home configuration
    qi = IK(Xi);
    dqi = [0; 0; 0];
    
    % Motion control with "Computed Torque Controller" (Task Space)
	taskSpaceMotionControl(pEqn,dpEqn,ddpEqn,HT,Jv,qi,dqi,tspan);
end

function [pEqn,dpEqn,ddpEqn] = taskSpaceTrajectory(Xi,Xf,dXi,dXf,tspan)
    syms q1 q2 q3 real;
    % Trajectory polynomial
    pEqn = sym(zeros(3,1)); 
    dpEqn = sym(zeros(3,1)); 
    ddpEqn = sym(zeros(3,1));
    % Create polynomial equation
    for n=1:3
        [pEqn(n),dpEqn(n),ddpEqn(n)] = generatePoly(Xi(n),Xf(n),dXi(n),dXf(n),tspan);
    end
end

% Generate polynomial for trajectory
function [qEqn,dqEqn,ddqEqn] = generatePoly(qi,qf,dqi,dqf,tspan)
    syms tt a0 a1 a2 a3 real
    % Polynomial trajectory
    qEqn = a0 + (a1*tt) + (a2*tt^2) + (a3*tt^3);
    dqEqn = diff(qEqn, tt);
    ddqEqn = diff(dqEqn, tt);
    % Plug in initial conditions
    qEqn_i = subs(qEqn, tt, tspan(1));
    dqEqn_i = subs(dqEqn, tt, tspan(1));
    qEqn_f = subs(qEqn, tt, tspan(2));
    dqEqn_f = subs(dqEqn, tt, tspan(2));
    % Solve for coefficients
    eqn = [qi==qEqn_i, dqi==dqEqn_i, qf==qEqn_f, dqf==dqEqn_f];
    [c0, c1, c2, c3] = solve(eqn, [a0,a1,a2,a3]);
    % Final trajectory equations
    qEqn = subs(qEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
    dqEqn = subs(dqEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
    ddqEqn = subs(ddqEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
end

function taskSpaceMotionControl(pEqn,dpEqn,ddpEqn,HT,Jv,qi,dqi,tspan)
    syms tt q1 q2 q3 real;
    global tme torque He
    
    tme = []; torque = []; He = [];
    dJv = diff(Jv);
    X0 = [qi; dqi];
    
    options = odeset('RelTol', 1e-4, 'AbsTol', [1e-4; 1e-4; 1e-4; 1e-4; 1e-4; 1e-4]);
	[T,X] = ode45(@(t,x)diffSolverTaskSpace(t,x,pEqn,dpEqn,ddpEqn,HT,Jv,dJv),tspan,X0,options);
    
    [eePos,eePosErr] = collectData(T,X,pEqn,HT);
    plotData(T,tme,X,torque,eePos,eePosErr,He);
    stickModel(HT,eePos,X);
end

function [eePos,eePosErr] = collectData(T,X,pEqn,HT)
    syms tt
    
    cnt = size(T,1);
    eePos = zeros(cnt,3);
    eePosErr = zeros(cnt,3);
    
    for n=1:cnt
        t = T(n);
        q = X(n,1:3)';
        % compute actual position
        p = FK(HT{end},q);
        % compute desired position
        pDes = eval(subs(pEqn,tt,t));
        % position error
        pErr = pDes - p;
        % update
        eePos(n,:) = p;
        eePosErr(n,:) = pErr;
    end
end

function dx = diffSolverTaskSpace(t,x,pEqn,dpEqn,ddpEqn,HT,Jv,dJv)
	global tme torque He  
    syms tt q1 q2 q3 real;
    
    xr = [0; 0; 50];
    K = 10*eye(3);
    Kp = [1000, 0, 0; 0, 1000, 0; 0, 0, 1000];
    Kv = [10, 0, 0; 0, 40, 0; 0, 0, 40];
    
    
    % Current Joint values
    curr_q = x(1:3,1);
	curr_dq = x(4:6,1);
    curr_Jv = eval(subs(Jv,[q1,q2,q3],curr_q'));
    curr_dJv = eval(subs(dJv,[q1,q2,q3],curr_q'));
    
    % Current end-effector values
    curr_p = FK(HT{end},curr_q);
    curr_dp = FVK(curr_Jv,curr_dq);
    
    % Compute current desired end-effector state vector:
    % [position, velocity, acceleration]
    des_p = eval(subs(pEqn, tt, t));
    des_dp = eval(subs(dpEqn, tt, t));
    des_ddp = eval(subs(ddpEqn, tt, t));
    
    % Operational Space Error
    pErr = des_p - curr_p;
    dpErr = des_dp - curr_dp;
    
    % Environment forces
    he = zeros(3,1);
    
    if ((curr_p(3)<50))
        he = inv((eye(3) + K*inv(Kp)))*K*(des_p-xr);
    end
    
    % Collect dynamic model
    [M,C,G] = dynamicModel(curr_q, curr_dq);
    
    % Compliance Control
    tau = G + transpose(curr_Jv)*Kp*pErr - transpose(curr_Jv)*Kv*curr_Jv*curr_dq;
    
    % Update state vector
    dx(1:3,1) = curr_dq;
    dx(4:6,1) = inv(M)*(tau - C.*curr_dq - G - transpose(curr_Jv)*he);
    
    torque = [torque; tau'];
    He = [He; he'];
    tme = [tme; t];
end


function [M,C,G] = dynamicModel(q, dq)
    m = [0.5, 0.5];
    a = 300; 
    b = 100; 
    c = 100;
    g = 9.81;
    M = [ m(2)*(b + q(2))^2,  0,  0; 0, m(2),  0; 0,  0, m(2)];
	C = [0; -(dq(1)^2*m(2)*(b + q(2))); 0];
    G = [0; 0; -g*m(2)];
end

function [T,linkPos] = hTrans()
    syms q1 q2 q3 real;
    q = [q1; q2; q3];
    a = 300; b = 100; c = 100;
    % Home configurations for each link
    M{1} = [1 0 0 0; 0 0 1 0; 0 -1 0 a; 0 0 0 1];
    M{2} = [1 0 0 0; 0 -1 0 b; 0 0 -1 a; 0 0 0 1];
    M{3} = [0 1 0 0; 1 0 0 b; 0 0 -1 a-c; 0 0 0 1];
    % Link velocities
    w1 = [0; 0; 1]; v1 = [0; 0; 0];
    w2 = [0; 0; 0]; v2 = [0; 1; 0];
    w3 = [0; 0; 0]; v3 = [0; 0; -1];
    % Concatenate velocites
    Ws = [w1 w2 w3]; 
    Vs = [v1 v2 v3];
    % Product of exponentials transformations
    ES = cell(3,1);
    % build transformations
    for i=1:3
        th = q(i);
        w = skew(Ws(:,i));
        R = eye(3)+(sin(th)*w)+((1-cos(th))*(w^2));
        v = (eye(3)*th+(1-cos(th))*w+(th-sin(th))*w^2)*Vs(:,i);
        ES{i} = [[R; 0, 0, 0], [v; 1]];
    end
    T = cell(3,1); linkPos = cell(3,1);
    T{1} = simplify(ES{1} * M{1}); 
    T{2} = simplify(ES{1} * ES{2} * M{2});
    T{3} = simplify(ES{1} * ES{2} * ES{3} * M{3});
    linkPos{1} = T{1}(1:3,4);
    linkPos{2} = T{2}(1:3,4);
    linkPos{3} = T{3}(1:3,4);
end

function W = skew(w)
    W = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

function J = Jacobian(linkPos)
    syms q1 q2 q3 real;
    J = simplify(jacobian(linkPos{end},[q1,q2,q3])); 
end

% Inverse kinematics
function [q] = IK(pos)
    q = zeros(3,1);
    x = pos(1); y = pos(2); z = pos(3);
    a = 300; 
    b = 100; 
    c = 100;
    q(3) = -z + a - c;
    q(1) = -atan2(x,y);
    q(2) = sqrt(x*x + y*y)- b;
end

% Compute joint velocities given end-effector velocity
function [dq] = IVK(J,dp)
    dq = pinv(J)*dp;
end

% Compute joint acceleration given end-effector acceleration
function [ddq] = IAK(J,dJ,dq,ddp)
    ddq = pinv(J)*(ddp - (dJ*dq));
end

% Forward kinematics
function [pos] = FK(T,q)
    syms q1 q2 q3 real;
    pos = eval(subs(T(1:3,4), [q1,q2,q3], q'));
end

% Compute end-effector velocity given joint velocities
function [dp] = FVK(J,dq)
    dp = J*dq;
end

% Compute end-effector acceleration given joint accelerations
function [ddp] = FAK(J,dJ,dq,ddq)
    ddp = dJ*dq + J*ddq;
end

function plotData(T,tme,X,torque,eePos,eePosErr,He)
    % Joint Position
    figure()
    for n=1:3
        plot(T, X(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('q_1', 'q_2', 'q_3')
    xlabel('time [sec]')
    ylabel('position [mm]')
    title(['Joint Position'])
    hold off
    
    % Joint Velocity
    figure()
    for n=4:6
        plot(T, X(:,n));
        hold on
    end
    xlim([0 T(end)])
    legend('dq_1', 'dq_2', 'dq_3')
    xlabel('time [sec]')
    ylabel('velocity [mm/sec]')
    title(['Joint Velocity'])
    hold off
    
    % Torque
	figure()
    for n=1:3
        plot(tme, torque(:,n));
        hold on
    end
    xlim([0 tme(end)])
    legend('q_1', 'q_2', 'q_3')
    xlabel('time [sec]')
    ylabel('Force [N*mm]')
    title(['Joint Torque'])
    hold off
    
    % End-Effector Position
    % x
    figure()
    subplot(1,3,1)
    plot(T,eePos(:,1))
    xlabel('time [sec]')
    ylabel('x [mm]')
    xlim([0 T(end)]);
    title('(x-axis)')
    %y
    subplot(1,3,2)
    plot(T,eePos(:,2))
    xlabel('time [sec]')
    ylabel('y [mm]')
    xlim([0 T(end)]);
    title('(y-axis)')
    %z
    subplot(1,3,3)
    plot(T,eePos(:,3))
    xlabel('time [sec]')
    ylabel('z [mm]')
    xlim([0 T(end)]);
    title('(z-axis)')
    sgtitle(['End-Effector Position'])
    
    % End-Effector Position Error
    % x
    figure()
    subplot(1,3,1)
    plot(T,eePosErr(:,1))
    xlabel('time [sec]')
    ylabel('x [mm]')
    xlim([0 T(end)]);
    title('(x-axis)')
    % y 
    subplot(1,3,2)
    plot(T,eePosErr(:,2))
    xlabel('time [sec]')
    ylabel('y [mm]')
    xlim([0 T(end)]);
    title('(y-axis)')
    % z
    subplot(1,3,3)
    plot(T,eePosErr(:,3))
    xlabel('time [sec]')
    ylabel('z [mm]')
    xlim([0 T(end)]);
    title('(z-axis)')
    sgtitle(['End-Effector Position Error'])
    drawnow
    
    % End-Effector Position Force
    % x
    figure()
    subplot(1,3,1)
    plot(tme,He(:,1))
    xlabel('time [sec]')
    ylabel('x [N*mm]')
    xlim([0 tme(end)]);
    title('(x-axis)')
    % y 
    subplot(1,3,2)
    plot(tme,He(:,2))
    xlabel('time [sec]')
    ylabel('y [N*mm]')
    xlim([0 tme(end)]);
    title('(y-axis)')
    % z
    subplot(1,3,3)
    plot(tme,He(:,3))
    xlabel('time [sec]')
    ylabel('z [N*mm]')
    xlim([0 tme(end)]);
    title('(z-axis)')
    sgtitle(['End-Effector Position Force'])
    drawnow
end

% Create stick model showing robot manipulator in its home position, and
% the trajectory of the end-effector
function [] = stickModel(HT,eePos,X)
    sz = size(eePos,1);
    traj = zeros(3,sz);
    q = X(1,1:3)';
    % Collect link positions in the home configuration
    p0 = [0; 0; 0];
    p1 = FK(HT{1},q);
    p2 = FK(HT{2},q);
    p3 = FK(HT{3},q);
    links = [p0, p1, p2, p3];
    % Plot the robot link
    figure()
    plot3(links(1,1:2), links(2,1:2), links(3,1:2));
    hold on
    plot3(links(1,2:3), links(2,2:3), links(3,2:3));
    hold on
    plot3(links(1,3:4), links(2,3:4), links(3,3:4));
    hold on
    % Collect the end-effector trajectory
    for n=1:sz
        traj(1:3,n) = eePos(n,:);
    end
    % Plot the end-effector trajectory
    plot3(traj(1,:), traj(2,:), traj(3,:),'--')
    legend('Link 1', 'Link 2', 'Link 3', 'End-Effector Trajectory')
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    title(['Robot Stick Model'])
    hold off
    drawnow    
end