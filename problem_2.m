function [] = problem_2()
    clear; close all; clc;
    syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
    
	q = [q1; q2; q3];
    dq = [dq1; dq2; dq3];
    ddq = [ddq1; ddq2; ddq3];

    T = hTrans(); % Forward Kinematics For Each Link
    Jv = Jacobian(T); % Jacobian For Each Link
    K = KE(Jv); 
    P = PE(T); % System Energies
    L = simplify(K - P); % Lagrange
    tau = [EL(L,q1,dq1); EL(L,q2,dq2); EL(L,q3,dq3)]; % Torque
    % Compact Dynamic Model
    M = mMatrix(tau);
    G = gMatrix(tau);
    C = cMatrix(tau,M,G);
    
    fprintf('\nInertia Matrix:\n');
    disp(M);
    fprintf('\nCoriolis and Centrifugal Matrix:\n');
    disp(C);
    fprintf('\n Gravity Matrix:\n');
    disp(G);

    Mc = eval(subs(M, [q1,q2,q3],[0,0,0]));
    Cc = eval(subs(C, [q1,q2,q3],[0,0,0])) ;
    Gc = eval(subs(G, [q1,q2,q3],[0,0,0]));
    
    fprintf('\nInertia Matrix at Home Position:\n');
    disp(Mc);
    fprintf('\nCoriolis and Centrifugal Matrix at Home Position:\n');
    disp(Cc);
    fprintf('\n Gravity Matrix at Home Position:\n');
    disp(Gc);
    
    el = simplify(expand(M*ddq + C.*dq + G));
    el = eval(subs(el, [q1,q2,q3],[0,0,0]));
    disp(['<strong>tau</strong> = ', char(el)]);
end

function T = hTrans()
    syms q1 q2 q3 real;
    L1 = 0.3; L2 = 0.3; L3 = 0.3;
    a = [0; L2; L3]; alpha = [sym(pi/2); 0; sym(-pi/2)];
    d = [L1; 0; 0]; theta = [q1; q2; q3];
    cnt = numel(a);
    T = cell(cnt,1); A = cell(cnt,1);
    T{1} = eye(4);
    for n=1:cnt
        A{n} = DH(alpha(n), a(n), d(n), theta(n));
    end
    T{1} = simplify(A{1});
    T{2} = simplify(A{1}*A{2});
    T{3} = simplify(A{1}*A{2}*A{3});
end

function Jv = Jacobian(T)
    syms q1 q2 q3 real;
    Jv = {};
    Jv{1} = jacobian(T{1}(1:3,4), [q1]);
    Jv{2} = jacobian(T{2}(1:3,4), [q1,q2]);
    Jv{3} = jacobian(T{3}(1:3,4), [q1,q2,q3]);
end

function K = KE(Jv)
    k = sym([]); m = 0.5*ones(3,1);
    V = velocity(Jv);
    k(1) = simplify((1/2)*m(1)*(dot(V{1}, V{1})));
    k(2) = simplify((1/2)*m(2)*(dot(V{2}, V{2})));
    k(3) = simplify((1/2)*m(3)*(dot(V{3}, V{3})));
    K = simplify(k(1) + k(2) + k(3));
end

function P = PE(T)
    p = sym([]); m = 0.5*ones(3,1);
    g = [0; 0; 9.8];
    p(1) = simplify(m(1)*g'*T{1}(1:3,4));
    p(2) = simplify(m(2)*g'*T{2}(1:3,4));
    p(3) = simplify(m(3)*g'*T{3}(1:3,4));
    P = simplify(p(1) + p(2) + p(3));
end

function V = velocity(Jv)
    syms dq1 dq2 dq3 real;
    V = {};
    V{1} = simplify(Jv{1}*dq1);
    V{2} = simplify(Jv{2}*[dq1; dq2]);
    V{3} = simplify(Jv{3}*[dq1; dq2; dq3]);
end

function tau = EL(L, q, dq)
    syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
    syms th1(t) th2(t) th3(t);

    L_ddq = diff(L, dq);
    L_t = subs(L_ddq,[q1 q2 q3 dq1 dq2 dq3], [th1, th2, th3, ...
        diff(th1(t),t), diff(th2(t),t), diff(th3(t), t)]);
    L_dt = diff(L_t, t);

    L1 = subs(L_dt,[th1, th2, th3, diff(th1(t),t), diff(th2(t),t), ...
     diff(th3(t),t), diff(th1(t),t,t), diff(th2(t),t,t), ...
     diff(th3(t),t,t)], [q1, q2, q3, dq1, dq2, dq3, ddq1, ddq2, ddq3]);
    L2 = diff(L, q);

    tau = simplify(L1 - L2);
end


function [A] = DH(alpha, a, d, theta)
    a11 = cos(theta);
    a12 = -1*sin(theta)*cos(alpha);
    a13 = sin(theta)*sin(alpha);
    a14 = a*cos(theta);
    a21 = sin(theta);
    a22 = cos(theta)*cos(alpha);
    a23 = -1*cos(theta)*sin(alpha);
    a24 = a*sin(theta);
    a32 = sin(alpha);
    a33 = cos(alpha);
    a34 = d;
    
    A = [a11, a12, a13, a14;
         a21, a22, a23, a24;
          0,  a32, a33, a34;
          0,   0,   0,   1];
end

function M = mMatrix(tau)
    M1 = mMatrixCol(tau(1));
    M2 = mMatrixCol(tau(2));
    M3 = mMatrixCol(tau(3));
    M = [M1; M2; M3];
    M = simplify(expand(M));
end

function M = mMatrixCol(tau)
    syms ddq1 ddq2 ddq3 real;
    ddq = [ddq1; ddq2; ddq3];
    m = sym([]);
    for n=1:3
        m(n) = simplify(tau - subs(tau,ddq(n),0))/ddq(n);
    end
    M = [m(1), m(2), m(3)];
end

function G = gMatrix(tau)
    syms dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
    dq = [dq1; dq2; dq3];
    ddq = [ddq1; ddq2; ddq3];

    %The gravity matrix is all terms not multiplied by dq or ddq.
    G = subs(tau, {ddq(1),ddq(2),ddq(3), dq(1),dq(2),dq(3)},{zeros(1,6)});
end

function C = cMatrix(tau,M,G)
    syms q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
    %The coriolis/cetripetal coupling vector is the result of subtracting
    %inertia and gravity portions from tau.
    C1 = simplify(expand(tau(1) - M(1,:)*[ddq1 ddq2 ddq3].' - G(1)));
    C2 = simplify(expand(tau(2) - M(2,:)*[ddq1 ddq2 ddq3].' - G(2)));
    C3 = simplify(expand(tau(3) - M(3,:)*[ddq1 ddq2 ddq3].' - G(3)));
    C = [C1;C2;C3];
end