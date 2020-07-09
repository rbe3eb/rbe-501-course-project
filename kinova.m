function [] = kinova()
    clc; clear; close all;
    
    % Homogeneous transform from DH parameters
    T = dh2Trans();
    [phi, theta, psi] = ZYZ(T);
    % compute jacobian
    J = computeJacobian(T);
    Ja = analyticalJacobian(J, phi, theta, psi);
    % Inverse kinematics
    desPos = [0, 0, 1];
    q = computeIK(T,J,desPos);

end

function [alpha,a,d,theta] = getDH()
    syms q1 q2 q3 q4 q5 q6 real
    % basic geometric parameters
    e2 = 0.0098;
    D1 = 0.2755; D2 = 0.4100; D3 = 0.2073;
    D4 = 0.7041; D5 = 0.7041; D6 = 0.1600;
    % DH parameters
    alpha = [sym(pi)/2; sym(pi); sym(pi)/2; sym(pi)/2; sym(pi)/2; sym(pi)];
    a = [0; D2; 0; 0; 0; 0];
    d = [D1; 0; e2; -(D3 + D4); 0; -(D5 + D6)];
    theta = [q1; q2; q3; q4; q5; q6];
end

function T = dh2Trans()
    % Get dh parameters
    [alpha,a,d,theta] = getDH();
    
    % Preallocate cells
    len = size(alpha,1);
    A = cell(1,len);
    T = cell(1,len);
    
    % Construct dh transformations
    for i=1:len
       A{i} = ...
       [cos(theta(i)) -sin(theta(i))*cos(alpha(i))  sin(theta(i))*sin(alpha(i)) a(i)*cos(theta(i))
        sin(theta(i))  cos(theta(i))*cos(alpha(i)) -cos(theta(i))*sin(alpha(i)) a(i)*sin(theta(i))
            0                 sin(alpha(i))               cos(alpha(i))               d(i)
            0                       0                            0                     1];
    end
    
    % Construct Forward kinematics
    T{1} = A{1};
    for i=2:len
        T{i} = simplify(T{i-1} * A{i});
    end
end

function J = computeJacobian(T)
    syms q1 q2 q3 q4 q5 q6
    Jv = simplify(jacobian(T{6}(1:3,4), [q1,q2,q3,q4,q5,q6]));
    Jw = [sym([0 0 1])', T{1}(1:3,3), T{2}(1:3,3), ...
        T{3}(1:3,3), T{4}(1:3,3), T{5}(1:3,3)];
    J = simplify([Jv; Jw]);
end

function Ja = analyticalJacobian(J, phi, theta, psi)
    B = [0 -sin(phi) sin(theta)*cos(phi)
         0  cos(phi) sin(theta)*sin(phi)
         1     0        cos(theta)];
     
    Ta = [  eye(3) zeros(3,3)
          zeros(3,3)  B];
    Ja = Ta\J;
end

function q = computeIK(T,J,desPos)
    syms q1 q2 q3 q4 q5 q6
    n = 0;
    desPos = [desPos 0 0 0];
    q = zeros(1,6);
    currPos = computeFK(T,zeros(1,6));
    currPos = [currPos 0 0 0];
    while n < 100 && sum(abs(desPos - currPos)) > 10^-10
        Jac = subsJ(J,q);
        posErr = pinv(Jac)*(desPos - currPos)';
        q = q + posErr';
        currPos = computeFK(T,q);
        currPos = [currPos 0 0 0];
        n = n + 1;
    end
    if sum(abs(desPos - currPos)) > 10^-10
        q = NaN*ones(1,6);
    end
end


function J = subsJ(J,q)
    syms q1 q2 q3 q4 q5 q6
    J = eval(subs(J, [q1,q2,q3,q4,q5,q6], q));
end

function pos = computeFK(T,q)
    syms q1 q2 q3 q4 q5 q6
    trans = eval(subs(T{6},[q1,q2,q3,q4,q5,q6], q));
    pos = trans(1:3,4)';
end

% Compute ZYZ-Euler angle transformation
function [phi, theta, psi] = ZYZ(T)
    phi = 0; theta = 0; psi = 0;
    % End-effector rotation matrix
    R = T{6}(1:3,1:3);
    if (R(1,3) ~= 0 && R(2,3) ~= 0)
        theta = atan2(R(3,3), sqrt(1-R(3,3)^2));
        phi = atan2(R(1,3), R(2,3));
        psi = atan2(-R(3,1), R(3,2));
    elseif (R(1,3) == 0 && R(2,3) == 0)
        if (R(3,3) == 1)
            theta = 0;
            phi = 0;
            psi = atan2(R(1,1), -R(1,2));
        elseif (R(3,3) == -1)
            theta = pi;
            phi = 0;
            psi = -atan2(-R(1,1), -R(1,2));
        end
    end
end

function eeVel = computeFVK(J,qVel)
    syms q1 q2 q3 q4 q5 q6
    J = subs(J, [q1,q2,q3,q4,q5,q6], qVel);
    eeVel = J*qVel;
end

function qVel = computeIVK(J,eeVel)
    syms q1 q2 q3 q4 q5 q6
    J = subs(J, [q1,q2,q3,q4,q5,q6], qVel);
    qVel = inv(J)*eeVel;
end