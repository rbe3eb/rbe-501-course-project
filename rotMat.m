function [R] = rotMat(T,q)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    R = eval(subs(T(1:3,1:3),[q1,q2,q3,q4,q5,q6,q7],q'));
end