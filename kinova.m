function [] = kinova()
    clc; clear;
    syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;
    % Compute homogeneous transformations
    T = hTran();
    % Compute Jacobian matrix at each link's tip
    [Jv, ~] = Jacobian(T);
    %[Ja] = analyticalJacobian(T{end},Jv{end},Jw{end});
    
    [M,C,G] = dynamicModel(T);
    %load('M.mat'); load('C.mat'); load('G.mat');
    % Motion Control: PD Plus Feed-Forward Controller
    Xi = [0; 200; 150] / 1000; 
    Xf = [200; 0; 200] / 1000;
    
    dXi = [0; 0; 0] / 1000; 
    dXf = [0; 0; 0] / 1000;
    
    ti = 0; tf = 10; % Initial and final time
    % Derive joint space trajectory polynomial
    [qEqn,dqEqn,ddqEqn] = jointSpaceTrajectory(Xi,Xf,dXi,dXf,ti,tf,T{end},Jv{end});
    jointSpaceMotionControl(qEqn,dqEqn,ddqEqn,T{end},Jv,M,C,G,[ti, tf]);
end

function [M,C,G] = dynamicModel(T)
    syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;

    % Link masses
    mass = [1.377, 1.1636, 1.1636, 0.930, 0.678, 0.678, 0.9 + 0.364];
    % Compute link centers of mass
    linkCOM = LinkCenterOfMass(T);
    % Collect inertia tensor
    I = inertialTensor();
    % Compute Jacobian matrix at each link' center of mass
    [Jv_COM, Jw_COM] = JacobianCOM(T, linkCOM);
    % Compute kinetic energy
	K = kineticEnergy(Jv_COM,Jw_COM,T,I,mass);
    % Compute potential energy
	P = potentialEnergy(linkCOM,mass);    
    L = simplify(K - P);
    % Compute torque
    tau = torque(L);
    % Compute Intertia matrix
    M = mMatrix(tau); 
    % Commpute Gravity matrix
    G = gMatrix(tau);
    % Centrifugal and coriolis matrix
    C = cMatrix(tau,M,G);
end

function [tau] = torque(L)
    syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;
    q = [q1; q2; q3; q4; q5; q6; q7];
    dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7];
    
    tau = sym(zeros(7,1));
    parfor n=1:7
        tau(n,1) = eulerLagrange(L,q(n),dq(n));
    end
end

function [] = jointSpaceMotionControl(qEqn,dqEqn,ddqEqn,HT,Jv,M,C,G,tspan) 
    global q dq eePos eePosErr torque;
    syms tt;
    q = []; dq = []; eePos = []; eePosErr = []; torque = [];
    
    % Initial position and velocity error
    X0 = ones(14,1);
    % Solve for joint errors
    [T,X] = ode45(@(t,x)diffSolverJointSpace(t,x,qEqn,dqEqn,ddqEqn,Jv,M,C,G,HT),tspan,X0);
    plotData(T, X, q, dq, toque, eePos);
end

function dx = diffSolverJointSpace(t,x,qEqn,dqEqn,ddqEqn,~,M,C,G,HT)
	global q dq eePos torque;   
    syms tt q1 q2 q3 q4 q5 q6 q7 real;
    
    Kp = 100*eye(7);
    Kv = 100*eye(7);
    dx = zeros(14,1);
    
    % Compute current desired motion -> [position, velocity, acceleration]
	qDes = eval(subs(qEqn, tt, t));
	dqDes = eval(subs(dqEqn, tt, t));
	ddqDes = eval(subs(ddqEqn, tt, t));
    
    % Obtain state variables -> [position err, velocity error]
	qErr = x(1:7,1);
	dqErr = x(8:14,1);
    
    curr_q = qDes - qErr;
    curr_dq = dqDes - dqErr;
    
    % Compute position and velocity
    curr_p = FK(HT,curr_q);
    
    % collect end-effector position
	q = [q, curr_q];
    dq = [dq, curr_dq];
    eePos = [eePos, curr_p];
    
    % Evaluate current dynamic model
    curr_M = eval(subs(M, [q1,q2,q3,q4,q5,q6,q7], curr_q'));
    curr_C = eval(subs(C, [q1,q2,q3,q4,q5,q6,q7], curr_q'));
    curr_G = eval(subs(G, [q1,q2,q3,q4,q5,q6,q7], curr_q'));
    % Evaluate desired dynamic model
    curr_Md = eval(subs(M, [q1,q2,q3,q4,q5,q6,q7], curr_q'));
    curr_Cd = eval(subs(C, [q1,q2,q3,q4,q5,q6,q7], curr_q'));
    curr_Gd = eval(subs(G, [q1,q2,q3,q4,q5,q6,q7], curr_q'));
    % Compute residual dynamics
    h = (curr_Md - curr_M)*ddqDes + (curr_Cd - curr_C)*dqDes + (curr_Gd - curr_G);
	tau = Kp*qErr + Kv*dqErr + curr_Md*ddqDes + curr_Cd*dqDes + curr_Gd;
    % Update state vector
    dx(1:7,1) = dqErr;
	dx(8:14,1) = inv(curr_M)*(-Kp*qErr - Kv*dqErr - curr_C*dqErr - h);
    torque = [torque, tau];
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
     
     % Frame 7 to interface
     T7i = [1, 0, 0, 0
         0, -1, 0, 0
         0, 0, -1, -0.0615
         0, 0, 0, 1];
     
     % Interface frame to gripper frame
     Tig = [0, 1, 0, 0
         1, 0, 0, 0
         0, 0, 1, 0.1493
         0, 0, 0, 1];
     
    % Frame 7 to gripper frame
    T6g = simplify(T67 * T7i * Tig);
    
    % Calculate homogeneous transformations
    T{1} = simplify(T01);
    T{2} = simplify(T01*T12);
    T{3} = simplify(T01*T12*T23);
    T{4} = simplify(T01*T12*T23*T34);
    T{5} = simplify(T01*T12*T23*T34*T45);
    T{6} = simplify(T01*T12*T23*T34*T45*T56);
    T{7} = simplify(T01*T12*T23*T34*T45*T56*T6g);
end

function [Jv_COM, Jw_COM] = JacobianCOM(T, linkCOM)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    % Linear velocity jacobian at each frame center of mass
    jv1 = simplify(jacobian(linkCOM(1:3,1), [q1,q2,q3,q4,q5,q6,q7]));
    jv2 = simplify(jacobian(linkCOM(1:3,2), [q1,q2,q3,q4,q5,q6,q7]));
    jv3 = simplify(jacobian(linkCOM(1:3,3), [q1,q2,q3,q4,q5,q6,q7]));
    jv4 = simplify(jacobian(linkCOM(1:3,4), [q1,q2,q3,q4,q5,q6,q7]));
    jv5 = simplify(jacobian(linkCOM(1:3,5), [q1,q2,q3,q4,q5,q6,q7]));
    jv6 = simplify(jacobian(linkCOM(1:3,6), [q1,q2,q3,q4,q5,q6,q7]));
    jv7 = simplify(jacobian(linkCOM(1:3,7), [q1,q2,q3,q4,q5,q6,q7]));
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
   Jv_COM = {jv1, jv2, jv3, jv4, jv5, jv6, jv7};
   Jw_COM = {jw1, jw2, jw3, jw4, jw5, jw6, jw7};
end

% Represent each link's jacobian at the tip of the link
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

function [Ja] = analyticalJacobian(T,Jv,Jw)
    J = [Jv, Jw];
    eul = rotm2eul(T(1:3,1:3), 'XYZ');
    r = eul(1);
    p = eul(2);
    y = eul(3);
    
    B = [1, 0, sin(p)
        0, cos(r), -cos(p)*sin(r)
        0, sin(r), cos(p)*cos(r)];
    
    Ta_inv = [eye(3), zeros(3,3)
        zeros(3,3), inv(B)];
    
    Ja = Ta_inv*J;
end

% Represent each link's center of mass wrt the base frame
function linkCOM = LinkCenterOfMass(T)
	linkCOM = sym(zeros(4,7));
    % Gripper center of mass wrt actuator 7 frame
    gripperCOM = [0; 0; -0.057-61.5];
    % Interface center of mass wrt actuator 7 frame
    interfaceCOM = [-0.000093; 0.000132; -0.022905];
    % Average center of mass of gripper and interface wrt actuator 7 frame
    endEffectorCOM = (gripperCOM + interfaceCOM) / 2;
    % Link of center of mass wrt the precedent joint reference frame
	L = sym([[-0.000648; -0.000166; 0.084487; 1], ...
        [-0.000023; -0.010364; -0.073360; 1], ...
        [-0.000044; -0.099580; -0.013278; 1], ...
        [-0.000044; -0.006641; -0.117892; 1], ...
        [-0.000018; -0.075478; -0.015006; 1], ...
        [0.000001; -0.009432; -0.063883; 1], ...
        [0.000001; -0.045483; -0.009650; 1], ...
        [endEffectorCOM; 1]]);
    % Link center of mass wrt the base frame
    linkCOM(1:4,1) = L(1:4,1);
    linkCOM(1:4,2) = simplify(T{1}*L(1:4,2));
    linkCOM(1:4,3) = simplify(T{2}*L(1:4,3));
    linkCOM(1:4,4) = simplify(T{3}*L(1:4,4));
    linkCOM(1:4,5) = simplify(T{4}*L(1:4,5));
    linkCOM(1:4,6) = simplify(T{5}*L(1:4,6));
    linkCOM(1:4,7) = simplify(T{6}*L(1:4,7));
    linkCOM = linkCOM(1:3,:);
end

function K = kineticEnergy(Jv,Jw,T,I,mass)
    syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
    
    dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7];
    m = cell(7,1);
    
    % Compute inertia matrix
    for n=1:7
        R = T{n}(1:3,1:3);
        a = simplify(mass(n)*Jv{n}'*Jv{n});
        b = simplify((Jw{n}'*R*I{n}*R'*Jw{n}));
    	m{n} = simplify(a + b);
    end
    M = simplify(m{1} + m{2} + m{3} + m{4} + m{5} + m{6} + m{7});
    % Compute kinetic energy
    K = simplify((1/2)*dq'*M*dq);
end

% Collect inertial tensors
function [I] = inertialTensor()
    I = cell(7,1);
    % Gripper moment of inertia wrt current frame
    IxxGripper = 0.00418; IxyGripper = 0; IxzGripper = 0;
    IyyGripper = 0.00518; IyzGripper = 0; IzzGripper = 0.00125;
    % Interface moment of inertia wrt current frame
    IxxInterface = 0.000214; IxyInterface = 0; IxzInterface = 0.000001;
    IyyInterface = 0.000223; IyzInterface = -0.000002; IzzInterface = 0.000240;
    % Combined moment of inertia of gripper and interface wrt current frame
    IxxEE = IxxGripper + IxxInterface; 
    IxyEE = IxyGripper + IxyInterface; 
    IxzEE = IxzGripper + IxzInterface;
    IyyEE = IyyGripper + IyyInterface;
    IyzEE = IyzGripper + IyzInterface;
    IzzEE = IzzGripper + IzzInterface;
    % Moment of inertia of each link wrt their center of mass
    Ixx = [0.004622, 0.004570 0.011088 0.010932 0.008147 0.001596 ...
        0.001641 IxxEE];
    Ixy = [0.000009, 0.000001 0.000005 0.000000 -0.000001 0.000000 ...
        0.000000 IxyEE];
    Ixz = [0.000060, 0.000002 0.000000 -0.000007 0.000000 0.000000 ...
        0.000000 IxzEE];
    Iyy = [0.004495, 0.004831 0.001072 0.011127 0.000631 0.001607 ...
        0.000410 IyyEE];
    Iyz = [0.000009, 0.000448 -0.000691 0.000606 -0.000500 0.000256 ...
        -0.000278 IyzEE];
    Izz = [0.002079, 0.001409 0.011255 0.001043 0.008316 0.000399 ...
        0.001641 IzzEE];
    Iyx = Ixy; 
    Izx = Ixz; 
    Izy = Iyz;
    % Create inertia tensors
    for n=1:7
        I{n} = [Ixx(n) Ixy(n) Ixz(n)
            Iyx(n) Iyy(n) Iyz(n)
            Izx(n) Izy(n) Izz(n)];
    end
end
    
function P = potentialEnergy(linkCOM,mass)
	g = [0; 0; 9.81];
    P = 0; p = sym(zeros(7,1));
    parfor n=1:7
        p(n) = simplify(mass(n)*g'*linkCOM(1:3,n));
    end
    P = simplify(p(1) + p(2) + p(3) + p(4) + p(5) + p(6) + p(7));
end

function M = mMatrix(tau)
    M = cell(7,1);
    parfor n=1:7
        M{n} = simplify(mMatrixCol(tau(n)));
    end
    M = [M{1}; M{2}; M{3}; M{4}; M{5}; M{6}; M{7}];
end

function M = mMatrixCol(taun)
    syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;
    ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6; ddq7];
    m = sym([]);
    parfor n=1:7
        m(n) = simplify(taun - subs(taun,ddq(n),0)/ddq(n));
    end
    M = [m(1), m(2), m(3), m(4), m(5), m(6), m(7)];
end

function G = gMatrix(tau)
    syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;
    dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7];
    ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6; ddq7];
    %The gravity matrix is all terms not multiplied by dq or ddq.
    G = simplify(subs(tau, {ddq(1),ddq(2),ddq(3),ddq(4),ddq(5),ddq(6),ddq(7), ...
        dq(1),dq(2),dq(3),dq(4),dq(5),dq(6),dq(7)},{zeros(1,14)}));
end

function c = cMatrix(tau,M,G)
    syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;
    %The coriolis/cetripetal coupling vector is the result of
    % subtracting inertia and gravity portions from tau.
    c = cell(7,1);
    
    parfor n=1:7
    	c{n} = simplify((tau(n) - M(n,:)*[ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ...
            ddq7].' - G(n)));
    end
    C = [c{1}; c{2}; c{3}; c{4}; c{5}; c{6}; c{7}];
end

function tau = eulerLagrange(L, q, dq)
    syms q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 dq7 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 ddq7 real;
    syms th1(t) th2(t) th3(t) th4(t) th5(t) th6(t) th7(t);
    
    L_ddq = diff(L, dq);
    L_t = subs(L_ddq,[q1 q2 q3 q4 q5 q6 q7 dq1 dq2 dq3 dq4 dq5 dq6 ...
        dq7], [th1, th2, th3, th4, th5, th6, th7, diff(th1(t),t), ...
        diff(th2(t),t), diff(th3(t),t), diff(th4(t),t), diff(th5(t),t), ...
        diff(th6(t), t), diff(th7(t), t)]);
    L_dt = diff(L_t, t);
    
    L1 = subs(L_dt,[th1, th2, th3, th4, th5, th6, th7, diff(th1(t),t), ...
        diff(th2(t),t), diff(th3(t),t), diff(th4(t),t), diff(th5(t),t), ...
        diff(th6(t),t), diff(th7(t),t), diff(th1(t),t,t), ...
        diff(th2(t),t,t), diff(th3(t),t,t), diff(th4(t),t,t), ...
        diff(th5(t),t,t), diff(th6(t),t,t), diff(th7(t),t,t),], ...
        [q1, q2, q3, q4, q5, q6, q7, dq1, dq2, dq3, dq4, dq5, dq6, dq7, ...
        ddq1, ddq2, ddq3, ddq4, ddq5, ddq6, ddq7]);
    
    L2 = diff(L, q);
    tau = simplify(L1 - L2);
end

% Generate trajectory in joint space
function [qEqn,dqEqn,ddqEqn] = jointSpaceTrajectory(Xi,Xf,dXi,dXf,ti,tf,T,Jv)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    % Trajectory polynomial
    qEqn = sym(zeros(7,1)); 
    dqEqn = sym(zeros(7,1)); 
    ddqEqn = sym(zeros(7,1));
    % initial and final joint positions
    qi = IK(Jv,T,Xi); %(Jv,T,pDes)
    qf = IK(Jv,T,Xf);
    % Jacobian
    curr_Jvi = eval(subs(Jv,[q1,q2,q3,q4,q5,q6,q7],qi'));
    curr_Jvf = eval(subs(Jv,[q1,q2,q3,q4,q5,q6,q7],qf'));
    % Initial and final joint velocities
    dqi = IVK(curr_Jvi,dXi); % (J,dp)
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

% Compute joint values given end-effector position
function curr_q = IK(Jv,T,pDes)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    q = [q1; q2; q3; q4; q5; q6; q7];
    curr_q = zeros(7,1);
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

function plotData(T, X, q, dq, torque, eePos)
    sz = size(T,1);
    % Plot Joint Position Error
    figure()
    for n=1:7
        plot(T, X(n,:));
        hold on
    end
    xlim([0 T(end)])
    legend('q_1', 'q_2', 'q_3', 'q_4', 'q_5', 'q_6', 'q_7')
    xlabel('time [sec]')
    ylabel('position [mm]')
    title('Joint Position Error')
    hold off
    
    % Plot Joint Velocity Error
    for n=8:14
        plot(T, X(n,:));
        hold on
    end
    xlim([0 T(end)])
    legend('dq_1', 'dq_2', 'dq_3', 'dq_4', 'dq_5', 'dq_6', 'dq_7')
    xlabel('time [sec]')
    ylabel('velocity [mm/sec]')
    title('Joint Velocity error')
    hold off
    
	figure()
    plot(T, torque(1, 1:sz));
    hold on
    plot(T, torque(2, 1:sz));
    hold on
    plot(T, torque(3, 1:sz));
    xlim([0 T(end)])
    legend('q_1', 'q_2', 'q_3')
    xlabel('time [sec]')
    ylabel('Force [N*mm]')
    title('Joint Torque (input)')
    hold off
    
    % Plot End-Effector Position 
    figure()
    subplot(1,3,1)
    plot(T,eePos(1:sz))
    xlabel('time [sec]')
    ylabel('x [mm]')
    xlim([0 T(end)]);
    title('End-Effector Position (x-axis)')
    
    subplot(1,3,2)
    plot(T,eePos(2:sz))
    xlabel('time [sec]')
    ylabel('y [mm]')
    xlim([0 T(end)]);
    title('End-Effector Position (y-axis)')
    
    subplot(1,3,3)
    plot(T,eePos(3:sz))
    xlabel('time [sec]')
    ylabel('z [mm]')
    xlim([0 T(end)]);
    title('End-Effector Position (z-axis)')
    sgtitle('End-Effector Position')    
end
    
