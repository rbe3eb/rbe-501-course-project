function [] = problem_1()
    clear all; close all; clc;
    syms g L1 L2 L3 I1 I2 I3 m1 m2 m3 q1 q2 q3 Lc1 Lc2 Lc3 dq1 dq2 dq3 dq1 dq2 dq3 ddq1 ddq2 ddq3 real;
    syms I1xx I2xx I3xx I1yy I2yy I3yy I1zz I2zz I3zz
    
    I1 = [I1xx I1yy I1zz];
    I2 = [I2xx I2yy I2zz];
    I3 = [I3xx I3yy I3zz];
    m = [m1 m2 m3];
    
    q = [q1; q2; q3];
    dq = [dq1; dq2; dq3];
    ddq = [ddq1; ddq2; ddq3];

    [T,R,Ri,com] = hTrans(); % Forward Kinematics For Each Link
    
    [~,Jvc,Jw] = Jacobian(T,com); % Jacobian For Each Link
    [~,D] = kineticEnergy(Jvc,Jw,R); 
    P = potentialEnergy(com);
    C = christoffel(D);
    G = gravity(P);
    
    % 1)
    fprintf('\n1. Inertia Matrix:\n');
    disp('D(q)= = ');
    disp(D);
    
    % 2)
    fprintf('\n2. Coriolis / Centrifugal Matrix:\n');
    disp('C(q,dq) = ');
    disp([C]);
    
    % 3)
	fprintf('\n3.1 Potential Energy:\n');
    disp('P = ');
    disp(P);
    fprintf('\n3.2 Gravity Matrix:\n');
    disp('G(q) = ');
    disp(G);
    
    % 4)
    el = simplify(expand(D*ddq + C*dq + G));
    fprintf('\n4. Compact Form: \n');
    disp(['tau = ', char(D), '*ddq + ', char(C), '*dq +', char(G), '\n']);

    % 5)
    ne = newtonEuler(T,R,Ri,com);
    fprintf('\n5. Newton-Euler Formulation:\n');
    disp('tau = ');
    disp(ne);
    
    % 6)
    el = simplify(expand(subs(el,[Lc1,Lc2,Lc3,I1,I2,I3],[L1,L2,L3,zeros(1,9)])),'Steps',50);
    el = simplify(expand(subs(el,[q',L1,L2,L3,m,g],[zeros(1,3),0.3,0.3,0.3,0.5*ones(1,3),9.8])),'Steps',50);
    fprintf('\n6. Numerical Dynamical Model:\n');
    disp(['tau = ',char(el)]);
    disp(newline);
    fprintf('*See that the result is identical to Problem 2-b from homework 3, which is located in the in the "Problem2_HW3" document\n')
end

function [G] = gravity(P)
    syms q1 q2 q3
    
    q = [q1; q2; q3];
    G = sym(zeros(3,1));
    for n=1:3
        G(n) = simplify(expand(diff(P,q(n))));
    end
end

function [C] = christoffel(D)
    syms q1 q2 q3 dq1 dq2 dq3 real
    
    q = [q1; q2; q3];
    dq = [dq1; dq2; dq3];
    C = sym(zeros(3,3));
    for k=1:3
        for j=1:3
            tempC = sym(zeros(3,1));
            for i=1:3
                a = simplify(diff(D(k,j),q(i)));
                b = simplify(diff(D(k,i),q(j)));
                c = simplify(diff(D(i,j),q(k)));
                tempC(i) = simplify(expand((1/2)*(a + b - c)))*dq(i);
            end
            C(k,j) = simplify(sum(tempC));
        end
    end
end

function [T,R,Ri,com] = hTrans()
    syms q1 q2 q3 L1 L2 L3 Lc1 Lc2 Lc3 real
    a = [0; L2; L3]; alpha = [sym(pi/2); 0; sym(-pi/2)];
    d = [L1; 0; 0]; theta = [q1; q2; q3];
    cnt = numel(a);
    T = cell(cnt,1); A = cell(cnt,1);
    R = cell(cnt,1); Ri= cell(cnt,1);
    com = cell(cnt,1);
    T{1} = eye(4);
    for n=1:cnt
        A{n} = DH(alpha(n), a(n), d(n), theta(n));
        Ri{n} = A{n}(1:3,1:3);
    end
    % Homogeneous transformations
    T{1} = simplify(A{1});
    T{2} = simplify(A{1}*A{2});
    T{3} = simplify(A{1}*A{2}*A{3});
    % Rotation matrices
    R{1} = T{1}(1:3,1:3);
    R{2} = T{2}(1:3,1:3);
    R{3} = T{3}(1:3,1:3);
    % Center of mass
    com{1} = subs(T{1}(1:4,4),L1,Lc1);
    com{2} = subs(T{2}(1:4,4),L2,Lc2);
    com{3} = subs(T{3}(1:4,4),L3,Lc3);
end

function [Jv,Jvc,Jw] = Jacobian(T,com)
    syms q1 q2 q3 real;
    
	Jv = {};
    Jvc = {};
    
    Jv{1} = [jacobian(T{1}(1:3,4), [q1]), zeros(3,2)];
    Jv{2} = [jacobian(T{2}(1:3,4), [q1,q2]), zeros(3,1)];
    Jv{3} = jacobian(T{3}(1:3,4), [q1,q2,q3]);
    
	Jvc{1} = [jacobian(com{1}(1:3), [q1]), zeros(3,2)];
    Jvc{2} = [jacobian(com{2}(1:3), [q1,q2]), zeros(3,1)];
    Jvc{3} = jacobian(com{3}(1:3), [q1,q2,q3]);
    
    Jw{1} = sym([[0; 0; 1], zeros(3,2)]);
    Jw{2} = [[0; 0; 1], T{1}(1:3,3), zeros(3,1)];
    Jw{3} = [[0; 0; 1], T{1}(1:3,3), T{2}(1:3,3)];
end

function [K,D] = kineticEnergy(Jvc,Jw,R)
    syms dq1 dq2 dq3 m1 m2 m3 real;
    syms I1xx I2xx I3xx I1yy I2yy I3yy I1zz I2zz I3zz
    
    I{1} = [I1xx, 0, 0; 0, I1yy, 0; 0, 0, I1zz];
    I{2} = [I2xx, 0, 0; 0, I2yy, 0; 0, 0, I2zz];
    I{3} = [I3xx, 0, 0; 0, I3yy, 0; 0, 0, I3zz];
    
    mass = [m1; m2; m3];
    dq = [dq1; dq2; dq3];
    d = cell(3,1);
    % Compute inertia matrix
    for n=1:3
        a = simplify(expand(mass(n)*transpose(Jvc{n})*Jvc{n}));
        b = simplify(expand((transpose(Jw{n})*R{n}*I{n}*transpose(R{n})*Jw{n})));
    	d{n} = simplify(expand(a + b));
    end
    D = simplify(expand(d{1} + d{2} + d{3}));
    % Compute kinetic energy
    K = simplify((1/2)*dq'*D*dq);
end


function P = potentialEnergy(com)
    syms dq1 dq2 dq3 I1 I2 I3 m1 m2 m3 g real;
    
    P = sym(0);
    mass = [m1; m2; m3];
    gi = [0; 0; g];
    for n=1:3
       P = simplify(P + mass(n)*gi'*com{n}(1:3)); 
    end
end


