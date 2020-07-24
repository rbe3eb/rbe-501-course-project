function [] = kinova()
    clc; clear;
    
    % Compute homogeneous transformations
    T = hTran();
    % Collect inertia tensor
    I = inertiaTensor();
    % Compute Jacobian matrix
    [Jv, Jw] = Jacobian(T);
	K = KE(Jv,Jw,T,I); % Kinetic Energy
	P = PE(T); % Potential Energy
    tau = torqueMatrix(K,P); % Torque Matrix
    M = mMatrix(tau); % Intertia matrix
    G = gMatrix(tau); % Gravity matrix
    C = cMatrix(tau,M,G); % Centrifugal and coriolis matrix
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
        0, 0, 1, 0
        -sin(q5), -cos(q5), 0, -0.1059
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
     
     % Frame 7 to gripper frame
     T7e = [1, 0, 0, 0
            0, 1, 0, 0
            0, 0, 1, 0.1250
            0, 0, 0, 1];
        
    % Calculate homogeneous transformations
    T{1} = simplify(T01);
    T{2} = simplify(T01*T12);
    T{3} = simplify(T01*T12*T23);
    T{4} = simplify(T01*T12*T23*T34);
    T{5} = simplify(T01*T12*T23*T34*T45);
    T{6} = simplify(T01*T12*T23*T34*T45*T56);
    T{7} = simplify(T01*T12*T23*T34*T45*T56*T67*T7e);
end

function [Jv, Jw] = Jacobian(T)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    linkCOM = LinkCenterOfMass(T);
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
   Jv = {jv1, jv2, jv3, jv4, jv5, jv6, jw7};
   Jw = {jw1, jw2, jw3, jw4, jw5, jw6, jv7};
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
	L = [[-0.000648; -0.000166; 0.084487], ...
        [-0.000023; -0.010364; -0.073360], ...
        [-0.000044; -0.099580; -0.013278], ...
        [-0.000044; -0.006641; -0.117892], ...
        [-0.000018; -0.075478; -0.015006], ...
        [0.000001; -0.009432; -0.063883], ...
        [0.000001; -0.045483; -0.009650], ...
        endEffectorCOM];
    % Link center of mass wrt the base frame
    linkCOM(1:4,1) = L(1:4,1);
    linkCOM(1:4,2) = simplify(T{2}*L(1:4,2));
    linkCOM(1:4,3) = simplify(T{3}*L(1:4,3));
    linkCOM(1:4,4) = simplify(T{4}*L(1:4,4));
    linkCOM(1:4,5) = simplify(T{5}*L(1:4,5));
    linkCOM(1:4,6) = simplify(T{6}*L(1:4,6));
    linkCOM(1:4,7) = simplify(T{7}*L(1:4,7));
    linkCOM = linkCOM(1:3,:);
end

function K = KE(Jv,Jw,T,I)
    syms dq1 dq2 dq3 dq4 dq5 dq6 dq7 real;
    
    dq = [dq1; dq2; dq3; dq4; dq5; dq6; dq7];
	M = sym(zeros(7,7));
    
    % Compute inertia matrix
    for n=1:7
    	M = M + (mass(n)*Jv{n}'*Jv{n}) + (Jw{n}'*T{n}(1:3,1:3)*I{n}* ...
            T{n}(1:3,1:3)'*Jw{n});
    end
    K = simplify((1/2)*dq'*M*dq);
end

% Collect inertial tensors
function I = inertialTensor()
    Ixx = [0.004622, 0.004570 0.046752 0.008292 0.001645 0.001685 ...
        0.000214];
    Ixy = [0.000009, 0.000001 -0.000009 -0.000001 0.000000 0.000000 ...
        0.000000];
    Ixz = [0.000060, 0.000002 0.000000 0.000000 0.000000 0.000000 ...
        0.000001];
    Iyy = [0.004495, 0.004831 0.000850 0.000628 0.001666 0.000400 ...
        0.000223];
    Iyz = [0.000009, 0.000448 -0.000098 0.000432 -0.000234 0.000255 ...
        0.000002];
    Izz = [0.002079, 0.001409 0.047188 0.008464 0.000389 0.001696 ...
        0.000240];
    Iyx = Ixy; Izx = Ixz; Izy = Iyz;
    I = cell(cnt,1);
    for n=1:cnt
        I{n} = [Ixx(n) Ixy(n) Ixz(n)
            Iyx(n) Iyy(n) Iyz(n)
            Izx(n) Izy(n) Izz(n)];
    end
end
    
function P = PE(T,mass,cnt)
	g = [0; 0; 9.81];
    r = LinkCenterOfMass(T,cnt);
    P = 0;
    for n=1:cnt
        P = P + mass(n)*g'*r(1:3,n);
    end
end

function tau =  torqueMatrix(K,P)
    syms q1 q2 q3 q4 q5 q6 ...
        dq1 dq2 dq3 dq4 dq5 dq6 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real;
    q = [q1; q2; q3; q4; q5; q6];
    dq = [dq1; dq2; dq3; dq4; dq5; dq6];
    ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6];
    
    L = simplify(K-P);
	dl_dq = derivative(L,q);
	dl_ddq = derivative(L,dq);
	dl_ddqdt = compute_velocity(dl_ddq,dq,ddq);
	tau = dl_ddqdt - dl_dq';
	tau = simplify(tau);
end

% derivative
function dldq = derivative(l,q)
    for i = 1:length(q)
        dldq(:,i) = diff(l,q(i));
    end
    dldq = simplify(dldq);
end

function vel = compute_velocity(p,q,dq)
    for i = 1:length(q)
        dpdq(:,i) = diff(p,q(i))*dq(i);
    end
    vel = simplify(expand(sum(dpdq,2)));
end

function M = mMatrix(tau)
    M1 = mMatrixCol(tau(1));
    M2 = mMatrixCol(tau(2));
    M3 = mMatrixCol(tau(3));
    M4 = mMatrixCol(tau(4));
    M5 = mMatrixCol(tau(5));
    M6 = mMatrixCol(tau(6));
    M = [M1; M2; M3; M4; M5; M6];
    M = simplify(expand(M));
end
    
function M = mMatrixCol(tau)
    syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real;
    ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6];
    m = cell(6,1);
    for n=1:6
        m{n} = simplify((tau - subs(tau,ddq(n),0))/ddq(n));
    end
    M = [m{1} m{2} m{3} m{4} m{5} m{6}];
end
    
function G = gMatrix(tau)
    syms dq1 dq2 dq3 dq4 dq5 dq6 ...
    ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real;
    dq = [dq1; dq2; dq3; dq4; dq5; dq6];
    ddq = [ddq1; ddq2; ddq3; ddq4; ddq5; ddq6];
        
    %The gravity matrix is all terms not multiplied by dq or ddq.
    G = subs(tau, {ddq(1),ddq(2),ddq(3),ddq(4),ddq(5),ddq(6), ...
        dq(1),dq(2),dq(3),dq(4),dq(5),dq(6)},{zeros(1,12)});   
end
    
function C = cMatrix(tau,M,G)
    syms q1 q2 q3 q4 q5 q6 ...
        dq1 dq2 dq3 dq4 dq5 dq6 ...
        ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real;
    %The coriolis/cetripetal coupling vector is the result of subtracting
    %inertia and gravity portions from tau.
    C1 = simplify(expand(tau(1) - M(1,:)*[ddq1 ddq2 ddq3 ddq4 ddq5 ...
        ddq6].' - G(1)));
    C2 = simplify(expand(tau(2) - M(2,:)*[ddq1 ddq2 ddq3 ddq4 ddq5 ...
        ddq6].' - G(2)));
    C3 = simplify(expand(tau(3) - M(3,:)*[ddq1 ddq2 ddq3 ddq4 ddq5 ...
        ddq6].' - G(3)));
    C4 = simplify(expand(tau(4) - M(4,:)*[ddq1 ddq2 ddq3 ddq4 ddq5 ...
        ddq6].' - G(4)));
    C5 = simplify(expand(tau(5) - M(5,:)*[ddq1 ddq2 ddq3 ddq4 ddq5 ...
        ddq6].' - G(5)));
    C6 = simplify(expand(tau(6) - M(6,:)*[ddq1 ddq2 ddq3 ddq4 ddq5 ...
        ddq6].' - G(6)));
    C = [C1;C2;C3;C4;C5;C6];
end
    
