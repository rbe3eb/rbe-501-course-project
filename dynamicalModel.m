function [M, C, G] = dynamicalModel()
    syms l1 l2 l3 l4 l5 l6
    syms q1 q2 q3 q4 q5 q6 real;
    syms dq1 dq2 dq3 dq4 dq5 dq6 real;
    syms ddq1 ddq2 ddq3 ddq4 ddq5 ddq6 real;
    syms m1 m2 m3 g real;
    
    % Mass of the link (grams)
    mass = [2292.668, 583.341, 1167.569, 593.781, 367.799, 367.799];
    
    % Length of the link (mm)
    length = [156.502, 159.758, 492.517, 248.258, 135.263, 135.263];
    
    % Coordinates of the center of mass (mm)
    rc = [[-0.301; 1.127; 70.634], [10.879; -0.01; 226.184], ...
        [-24.195; -0.002; 480.492], [16.852; 0.002; 760.264], ...
        [-4.545; 0.015; 956.344], [24.045; 0.015; 1036.656]];
    
    % Moment of Inertia at Center of Mass (g*mm^2)
	Ixx = [5.263E+06, 1.902E+06, 2.870E+07, 3.912E+06, 7.368E+05, 7.368E+05];
	Ixy = [2.418E+04, 9.402, -7.653, 17.362, -8.345,  8.345];
	Ixz = [2.709E+04, -2.331E+05, -43.188, 1.489E+05, 1.510E+05, 1.510E+05];
	Iyx = [2.418E+04, 9.402, -7.653, 17.362, -8.345, 8.345];
	Iyy = [5.131E+06, 1.778E+06, 2.819E+07, 3.727E+06, 7.492E+05, 7.492E+05]; 
	Iyz = [-2.362E+05, -98.997, -867.047, 334.219, 63.799, -63.799];
	Izx = [2.709E+04, -2.331E+05, -43.188, 1.489E+05, 1.510E+05, 1.510E+05]; 
	Izy = [-2.362E+05, -98.997, -867.047, 334.219, 63.799, -63.799];
	Izz = [2.594E+06, 7.121E+05, 8.923E+05, 4.649E+05, 3.035E+05, 3.035E+05];
    
    
    
    %
    a = [0 l2 l3];
    alpha = [sym(pi)/2 0 -sym(pi)/2];
    d = [l1 0 0];
    theta = [q1 q2 q3];

    T01 = DH(alpha(1), a(1), d(1), theta(1));
    T12 = DH(alpha(2), a(2), d(2), theta(2));
    T23 = DH(alpha(3), a(3), d(3), theta(3));

    T02 = simplify(T01 * T12);
    T03 = simplify(T02 * T23);

    % plot3([0 0],[0 0],[-100 0],'-bo','LineWidth',10,...
    % 'MarkerEdgeColor','k',...
    % 'MarkerFaceColor','g',...
    % 'MarkerSize',1);
    % 
    % grid on;
    % hold on;


    t = [0 0 0];
    T01_temp = subs(T01, [l1 l2 l3 q1 q2 q3], [1 1 1 t(1) t(2) t(3)]);
    T02_temp = subs(T02, [l1 l2 l3 q1 q2 q3], [1 1 1 t(1) t(2) t(3)]);
    T03_temp = subs(T03, [l1 l2 l3 q1 q2 q3], [1 1 1 t(1) t(2) t(3)]);

    x = [0 T01_temp(1,4) T02_temp(1,4) T03_temp(1,4)];
    y = [0 T01_temp(2,4) T02_temp(2,4) T03_temp(2,4)];
    z = [0 T01_temp(3,4) T02_temp(3,4) T03_temp(3,4)];

    % plot3(x,y,z,'-ro','LineWidth',4,...
    % 'MarkerEdgeColor','k',...
    % 'MarkerFaceColor','g',...
    % 'MarkerSize',10);

    J_l1 = jacobian(T01(1:3,4),[q1]);
    J_l2 = jacobian(T02(1:3,4),[q1 q2]);
    J_l3 = jacobian(T03(1:3,4),[q1 q2 q3]);

    J_l1 = simplify(J_l1)
    J_l2 = simplify(J_l2)
    J_l3 = simplify(J_l3)
    return

    v_m1 = J_l1 * dq1;
    v_m2 = J_l2 * [dq1 ; dq2];
    v_m3 = J_l3 * [dq1 ; dq2 ; dq3];
    v_m1 = simplify(v_m1);
    v_m2 = simplify(v_m2);
    v_m3 = simplify(v_m3);

    K1 = 0.5 * m1 * (v_m1.' * v_m1);
    K2 = 0.5 * m2 * (v_m2.' * v_m2);
    K3 = 0.5 * m3 * (v_m3.' * v_m3);
    K1 = simplify(K1);
    K2 = simplify(K2);
    K3 = simplify(K3);

    g = 9.8;

    P1 = m1 * g * T01(3,4);
    P2 = m2 * g * T02(3,4);
    P3 = m3 * g * T03(3,4);

    P1 = simplify(P1);
    P2 = simplify(P2);

    P3 = simplify(P3);

    K = K1 + K2 + K3;
    P = P1 + P2 + P3;


    L = simplify( K - P );
    q = [q1 q2 q3];
    dq = [dq1 dq2 dq3];
    tau = cell(3,1);
    for n=1:3
        tau{n} = lagrange(L,q(n),dq(n));
    end
    syms th1(t) th2(t) th3(t);
    
    A1 = diff(L,dq1);
    A1t = subs(A1,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
    dA1t = diff(A1t, t);
    A1 = subs(dA1t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t) ...
        diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
    B1 = diff(L,q1);

    A2 = diff(L,dq2);
    A2t = subs(A2,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
    dA2t = diff(A2t, t);
    A2 = subs(dA2t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t),t) ...
        diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);
    B2 = diff(L,q2);

    A3 = diff(L,dq3);
    A3t = subs(A3,[q1 q2 q3 dq1 dq2 dq3], [th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t)]);
    dA3t = diff(A3t, t);
    A3 = subs(dA3t,[th1 th2 th3 diff(th1(t),t) diff(th2(t),t) diff(th3(t), t) ...
        diff(th1(t),t,t) diff(th2(t),t,t) diff(th3(t), t,t)],[q1 q2 q3 dq1 dq2 dq3 ddq1 ddq2 ddq3]);

    B3 = diff(L,q3);

    Tau_l1 = A1 - B1;
    Tau_l2 = A2 - B2;
    Tau_l3 = A3 - B3;

    M11 = simplify(Tau_l1 - subs(Tau_l1,ddq1,0)) /ddq1;
    M12 = simplify(Tau_l1 - subs(Tau_l1,ddq2,0)) /ddq2;
    M13 = simplify(Tau_l1 - subs(Tau_l1,ddq3,0)) /ddq3;

    M21 = simplify(Tau_l2 - subs(Tau_l2,ddq1,0)) /ddq1;
    M22 = simplify(Tau_l2 - subs(Tau_l2,ddq2,0)) /ddq2;
    M23 = simplify(Tau_l2 - subs(Tau_l2,ddq3,0)) /ddq3;

    M31 = simplify(Tau_l3 - subs(Tau_l3,ddq1,0)) /ddq1;
    M32 = simplify(Tau_l3 - subs(Tau_l3,ddq2,0)) /ddq2;
    M33 = simplify(Tau_l3 - subs(Tau_l3,ddq3,0)) /ddq3;

    M = [M11 M12 M13;
        M21 M22 M23;
        M31 M32 M33];
    M = simplify(M);

    G11 = subs(Tau_l1, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
    G21 = subs(Tau_l2, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
    G31 = subs(Tau_l3, [dq1 dq2 dq3 ddq1 ddq2 ddq3], [0 0 0 0 0 0]);
    G = simplify([G11; G21; G31]);

    C11 = Tau_l1 - (M(1,:) * [ddq1 ddq2 ddq3].' + G11);
    C21 = Tau_l2 - (M(2,:) * [ddq1 ddq2 ddq3].' + G21);
    C31 = Tau_l3 - (M(3,:) * [ddq1 ddq2 ddq3].' + G31);
    C = simplify([C11; C21; C31]);
    
    % Physical properties
    l = [0.3,0.3,0.3];
    m = [0.5,0.5,0.5];
    % Substitute physical properties
    M = subs(M, [l1,l2,l3,m1,m2,m3], [m,l]);
    C = subs(C, [l1,l2,l3,m1,m2,m3], [m,l]);
    G = subs(G, [l1,l2,l3,m1,m2,m3], [m,l]);
    % Create struct containing dynamical model
    DM = struct;
    DM.M = M;
    DM.C = C;
    DM.G = G;
    
    T01 = subs(T01, [l1,l2,l3],l);
    T02 = subs(T02, [l1,l2,l3],l);
    T03 = subs(T03, [l1,l2,l3],l);
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