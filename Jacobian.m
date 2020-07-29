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