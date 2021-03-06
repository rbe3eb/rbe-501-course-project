function [] = control()
    clc; clear;
    syms q1 q2 q3 q4 q5 q6 q7 real;
    
    % Compute homogeneous transformations
    HT = hTran(); [Jv, ~] = Jacobian(HT);
    
    % Initial and final time
    ti = 0; tf = 4;
    
    % Initial joint values and velocites
    qi = [0; 0; 0; 0; 0; 0; 0];
    dqi = [0; 0; 0; 0; 0; 0; 0];
    
    % Final joint positions and velocities
    Xf = {[-0.25; 0.55; 0.44], [-0.05; 0.55; 0.44], [-0.05; 0.55; 0.62], ...
        [0.3; 0.5; 0.26]};
    dXf = {zeros(3,1), zeros(3,1), zeros(3,1), zeros(3,1)};
    
    % Derive joint space trajectory polynomial
    [qEqn,dqEqn,ddqEqn] = jointSpaceTrajectory(HT{end},Jv{end},qi,dqi,Xf{1},dXf{1},ti,tf);
    
    % Simulate robot motion
    %jointSpaceMotionControl(qEqn,dqEqn,ddqEqn,HT{end},Jv{end},M,C,G,[ti, tf]);
end

% Generate trajectory in joint space
function [qEqn,dqEqn,ddqEqn] = jointSpaceTrajectory(HT,Jv,qi,dqi,Xf,dXf,ti,tf)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    
    % Trajectory polynomial
    qEqn = sym(zeros(7,1)); 
    dqEqn = sym(zeros(7,1)); 
    ddqEqn = sym(zeros(7,1));
    
    % final joint positions
    qf = IK(HT,Jv,Xf,qi); % (T,Jv,pDes,curr_q)
    
    % Evaluate the Jacobian
    curr_Jvf = eval(subs(Jv,[q1,q2,q3,q4,q5,q6,q7],qf'));
    
    % Initial and final joint velocities
    dqf = IVK(curr_Jvf,dXf);
    
    % Create polynomial equations for each joint
    for n=1:7
        [qEqn(n),dqEqn(n),ddqEqn(n)] = generatePoly(qi(n),qf(n),dqi(n),dqf(n),ti,tf);
    end
end


% Generate polynomial for trajectory
function [qEqn,dqEqn,ddqEqn] = generatePoly(qi,qf,dqi,dqf,ti,tf)
    syms tt a0 a1 a2 a3 real
    % Polynomial trajectory
    qEqn = a0 + (a1*tt) + (a2*tt^2) + (a3*tt^3);
    dqEqn = diff(qEqn, tt);
    ddqEqn = diff(dqEqn, tt);
    % Plug in initial conditions
    qEqn_i = subs(qEqn, tt, ti);
    dqEqn_i = subs(dqEqn, tt, ti);
    qEqn_f = subs(qEqn, tt, tf);
    dqEqn_f = subs(dqEqn, tt, tf);
    % Solve for coefficients
    eqn = [qi==qEqn_i, dqi==dqEqn_i, qf==qEqn_f, dqf==dqEqn_f];
    [c0, c1, c2, c3] = solve(eqn, [a0,a1,a2,a3]);
    % Final trajectory equations
    qEqn = subs(qEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
    dqEqn = subs(dqEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
    ddqEqn = subs(ddqEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
end

function [] = jointSpaceMotionControl(qEqn,dqEqn,ddqEqn,HT,Jv,qi,dqi,ti,tf) 
    global eePos eePosErr torque tme;
    eePos = []; eePosErr = []; torque = []; tme = [];
    
    % Initial position and velocity error
    X0 = [qi; dqi];
    
    % Solve for joint position and velocity errors
    options = odeset('RelTol', 1e-4, 'AbsTol', [1e-4; 1e-4; 1e-4; 1e-4; 1e-4; 1e-4]);
    [T,X] = ode45(@(t,x)diffSolverJointSpace(t,x,qEqn,dqEqn,ddqEqn,HT,Jv,method),[ti tf],X0,options);
end

function dx = diffSolverJointSpace(t,x,qEqn,dqEqn,ddqEqn,HT,~)
	global eePos eePosErr torque tme;   
    syms tt;
    
    dx = zeros(14,1);
    Kp = 100*eye(7);
    Kv = 100*eye(7);
    
    % Compute current desired motion -> [position, velocity, acceleration]
    qDes = eval(subs(qEqn{1}, tt, t));
    dqDes = eval(subs(dqEqn{1}, tt, t));
    ddqDes = eval(subs(ddqEqn{1}, tt, t));
    
    % Obtain state variables -> [joint angle error, joint velocity error]
	curr_q = x(1:7,1);
	curr_dq = x(8:14,1);
    
    % Compute error in the position and velocity
    qErr = qDes - curr_q;
    dqErr = dqDes - curr_dq;
    
    % Compute end-effector position and end-effector position error
    curr_eePos = FK(HT{end}, curr_q);
    pDes = FK(HT{end}, qDes);
    curr_eePosErr = pDes - curr_eePos;
    
    % collect end-effector position and position error
    eePos = [eePos; curr_eePos'];
    eePosErr = [eePosErr; curr_eePosErr'];
    tme = [tme; t];
    
    % Collect dynamic model
    [M,C,G] = dynamicModel(curr_q, curr_dq);
    
    tau = M*(ddqDes + Kv*dqErr + Kp*qErr) + C.*curr_dq + G;

    % Update state vector
    dx(1:7,1) = curr_dq;
    dx(8:14,1) = inv(M)*(tau - C.*curr_dq - G);
    torque = [torque; tau'];
end

function T = hTran()
	syms q1 q2 q3 q4 q5 q6 q7 real
    T = cell(7,1);
    
    % Base to frame 1
    T01 = [cos(q1), -sin(q1), 0, 0
           -sin(q1), -cos(q1), 0, 0
           0, 0, -1, 0.1564
           0, 0, 0, 1];
       
    % Frame 1 to frame 2
	T12 = [cos(q2), -sin(q2), 0, 0
        0, 0, -1, 0.0054
        sin(q2), cos(q2), 0, -0.1284
        0, 0, 0, 1];
    
    % Frame 2 to frame 3
	T23 = [cos(q3), -sin(q3), 0, 0
        0, 0, 1, -0.2104
        -sin(q3), -cos(q3), 0, -0.0064
        0, 0, 0, 1];
    
    % Frame 3 to frame 4
	T34 = [cos(q4), -sin(q4), 0, 0
        0, 0, -1, 0.0064
        sin(q4), cos(q4), 0, -0.2104
        0, 0, 0, 1];

    % Frame 4 to frame 5
	T45 = [cos(q5), -sin(q5), 0, 0
        0, 0, 1, -0.2084
        -sin(q5), -cos(q5), 0, -0.0064
        0, 0, 0, 1];
    
    % Frame 5 to frame 6
     T56 = [cos(q6), -sin(q6), 0, 0
         0, 0, -1, 0
         sin(q6), cos(q6), 0, -0.1059
         0, 0, 0, 1];
     
     % Frame 6 to frame 7
     T67 = [cos(q7), -sin(q7), 0, 0
         0, 0, 1, -0.1059
         -sin(q7), -cos(q7), 0, 0
         0, 0, 0, 1];
     
    T7_ee = [1, 0, 0, 0
        0, -1, 0, 0
        0, 0, -1, -0.0615
        0, 0, 0, 1];

    Tee_g = [1, 0, 0, 0
        0, 1, 0, 0
        0, 0, 1, 0.1250
        0, 0, 0, 1];

    % Frame 7 to gripper frame
    T6g = simplify(T67 * T7_ee * Tee_g);
    
    % Calculate homogeneous transformations
    T{1} = simplify(T01);
    T{2} = simplify(T01*T12);
    T{3} = simplify(T01*T12*T23);
    T{4} = simplify(T01*T12*T23*T34);
    T{5} = simplify(T01*T12*T23*T34*T45);
    T{6} = simplify(T01*T12*T23*T34*T45*T56);
    T{7} = simplify(T01*T12*T23*T34*T45*T56*T6g);
end

function [Jv, Jw] = Jacobian(T)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    % Linear velocity jacobian at each frame center of mass
    jv1 = simplify(jacobian(T{1}(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    jv2 = simplify(jacobian(T{2}(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    jv3 = simplify(jacobian(T{3}(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    jv4 = simplify(jacobian(T{4}(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    jv5 = simplify(jacobian(T{5}(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    jv6 = simplify(jacobian(T{6}(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    jv7 = simplify(jacobian(T{7}(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    % Angular velocity jacobian at each center of mass
    jw1 = [[0; 0; 1], zeros(3,6)];
    jw2 = [[0;0;1], T{1}(1:3,3), zeros(3,5)];
    jw3 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), zeros(3,4)];
    jw4 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), T{3}(1:3,3), zeros(3,3)];
    jw5 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), T{3}(1:3,3), ...
        T{4}(1:3,3), zeros(3,2)];
    jw6 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), T{3}(1:3,3), ...
        T{4}(1:3,3), T{5}(1:3,3), zeros(3,1)];
    jw7 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), T{3}(1:3,3), ...
        T{4}(1:3,3), T{5}(1:3,3), T{6}(1:3,3)];
   % Combine jacobian matrices
   Jv = {jv1, jv2, jv3, jv4, jv5, jv6, jv7};
   Jw = {jw1, jw2, jw3, jw4, jw5, jw6, jw7};
end


% Compute joint values given end-effector position
function curr_q = IK(T,Jv,pDes,curr_q)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    q = [q1; q2; q3; q4; q5; q6; q7];

    curr_p = FK(T,curr_q);
	n = 1;
    epsilon = 0.00001;
    while ((norm(pDes - curr_p)) > epsilon && n < 400)
        if isequal(mod(n,20), 0)
            curr_q = 0 + (pi-0).*rand(7,1);
            curr_p = FK(T,curr_q);
        end
        % Evaluate jacobian given current joint values
        curr_Jv = eval(subs(Jv, q, curr_q));
        % Determine desired change in joint values
        qDelta = pinv(curr_Jv)*(pDes-curr_p);
        % Update joint values
        curr_q = curr_q + qDelta;
        % Evaluate position given current joint values
        curr_p = FK(T,curr_q);
        n = n + 1;
    end
end

% Compute joint velocities given end-effector velocity
function [dq] = IVK(J,dp)
    dq = pinv(J)*dp;
end

% Compute joint acceleration given end-effector acceleration
function [ddq] = IAK(J,dJ,dq,ddp)
    ddq = pinv(J)*(ddp - (dJ*dq));
end

% Compute end-effector position given joint values
function [p] = FK(T,q)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    p = eval(subs(T(1:3,4),[q1,q2,q3,q4,q5,q6,q7],q'));
    p = p(1:3,1);
end

% Compute end-effector velocity given joint velocities
function [dp] = FVK(J,dq)
    dp = J*dq;
end

% Compute end-effector acceleration given joint accelerations
function [ddp] = FAK(J,dJ,dq,ddq)
    ddp = dJ*dq + J*ddq;
end
    
