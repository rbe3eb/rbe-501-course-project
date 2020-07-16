function [] = kinova()
    clc; clear;
    
    mass = [1.697, 1.377, 1.262, 0.930, 0.678, 0.678, 0.364];
    cnt = 6;
    
    T = transformationMatrix(cnt); % Transformation matrix
    [Jv, Jw] = Jacobian(T,cnt); % Jacobian matrix
	K = KE(Jv,Jw,T,mass,cnt); % Kinetic Energy
	P = PE(T,mass,cnt); % Potential Energy
    tau = torqueMatrix(K,P); % Torque Matrix
    M = mMatrix(tau); % Intertia matrix
    G = gMatrix(tau); % Gravity matrix
    C = cMatrix(tau,M,G); % Centrifugal and coriolis matrix
end

function T = transformationMatrix(cnt)
	syms q1 q2 q3 q4 q5 q6 real
    T = cell(7,1);

    T01 = [cos(q1), -sin(q1), 0, 0
           -sin(q1), -cos(q1), 0, 0
           0, 0, -1, 0.1564
           0, 0, 0, 1];

     T12 = [cos(q2), -sin(q2), 0, 0
             0, 0, -1, 0.0054
             sin(q2), cos(q2), 0, -0.1284
             0, 0, 0, 1];

     T23 = [cos(q3), -sin(q3), 0, 0
             -sin(q3), -cos(q3), 0, -0.410
             0, 0, -1, 0
             0, 0, 0, 1];

     T34 = [cos(q4), -sin(q4), 0, 0
             0, 0, -1, 0.2084
             sin(q4), cos(q4), 0, -0.0064
             0, 0, 0, 1];

     T45 = [cos(q5), -sin(q5), 0, 0
             0, 0, 1, 0
             -sin(q5), -cos(q5), 0, -0.1059
             0, 0, 0, 1];
     T56 = [cos(q6), -sin(q6), 0, 0
            0, 0, -1, 0.1059
            sin(q6), cos(q6), 0, 0
            0, 0, 0, 1];
     T67 = [-1, 0, 0, 0
            0, 1, 0, 0
            0, 0, -1, -0.615
            0, 0, 0, 1];
        
    T{1} = simplify(T01);
    T{2} = simplify(T01*T12);
    T{3} = simplify(T01*T12*T23);
    T{4} = simplify(T01*T12*T23*T34);
    T{5} = simplify(T01*T12*T23*T34*T45);
    T{6} = simplify(T01*T12*T23*T34*T45*T56);
    T{7} = simplify(T01*T12*T23*T34*T45*T56*T67);
end

function [Jv, Jw] = Jacobian(T,cnt)
    syms q1 q2 q3 q4 q5 q6 real;
    r = COM(T,cnt);
    jv1 = simplify(jacobian(r(1:3,1), [q1,q2,q3,q4,q5,q6]));
    jv2 = simplify(jacobian(r(1:3,2), [q1,q2,q3,q4,q5,q6]));
    jv3 = simplify(jacobian(r(1:3,3), [q1,q2,q3,q4,q5,q6]));
    jv4 = simplify(jacobian(r(1:3,4), [q1,q2,q3,q4,q5,q6]));
    jv5 = simplify(jacobian(r(1:3,5), [q1,q2,q3,q4,q5,q6]));
    jv6 = simplify(jacobian(r(1:3,6), [q1,q2,q3,q4,q5,q6]));

    jw1 = [[0; 0; 1], zeros(3,5)];
    jw2 = [[0;0;1], T{1}(1:3,3), zeros(3,4)];
    jw3 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), zeros(3,3)];
    jw4 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), T{3}(1:3,3), zeros(3,2)];
    jw5 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), T{3}(1:3,3), ...
        T{4}(1:3,3), zeros(3,1)];
    jw6 = [[0;0;1], T{1}(1:3,3), T{2}(1:3,3), T{3}(1:3,3), ...
        T{4}(1:3,3), T{5}(1:3,3)];
      
   Jv = {jv1, jv2, jv3, jv4, jv5, jv6};
   Jw = {jw1, jw2, jw3, jw4, jw5, jw6};
end

function r = COM(T,cnt)
	r = sym(zeros(4,cnt));
	L = [[-0.000648; -0.000166; 0.084487; 1], ...
        [-0.000023; -0.010364; -0.073360; 1], ...
        [0.000035; -0.208207; -0.018890; 1], ...
        [0.000018; 0.076168; -0.013970; 1], ...
        [-0.000001; 0.008466; -0.062937; 1], ...
        [-0.000001; 0.046429; -0.008704; 1], ...
        [-0.000093; 0.000132; -0.022905; 1]];
    r(1:4,1) = simplify(T{1}*L(1:4,1));
    r(1:4,2) = simplify(T{2}*L(1:4,2));
    r(1:4,3) = simplify(T{3}*L(1:4,3));
    r(1:4,4) = simplify(T{4}*L(1:4,4));
    r(1:4,5) = simplify(T{5}*L(1:4,5));
    r(1:4,6) = simplify(T{6}*L(1:4,6));
    r(1:4,7) = simplify(T{7}*L(1:4,7));
    r = r(1:3,:);
end

function K = KE(Jv,Jw,T,mass,cnt)
    syms q1 q2 q3 q4 q5 q6 ...
        dq1 dq2 dq3 dq4 dq5 dq6 real;
    I = inertialParams(cnt);
    dq = [dq1; dq2; dq3; dq4; dq5; dq6];
    M = inertiaMatrix(Jv,Jw,I,T,mass,cnt);
    K = simplify((1/2)*dq'*M*dq);
end

function I = inertialParams(cnt)
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

function M = inertiaMatrix(Jv,Jw,I,T,mass,cnt)
	M = sym(zeros(cnt,cnt));
    for n=1:cnt
    	M = M + (mass(n)*Jv{n}'*Jv{n}) + (Jw{n}'*T{n}(1:3,1:3)*I{n}* ...
            T{n}(1:3,1:3)'*Jw{n});
        %M = simplify(M);
    end
end
    
function P = PE(T,mass,cnt)
	g = [0; 0; 9.81];
    r = COM(T,cnt);
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
    
