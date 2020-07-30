function [gripper] = gripperPos(Tg,eePos,q)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    
    cnt = size(eePos,1);
    gripper = (zeros(cnt,3));
    Tg = eval(subs(Tg,[q1,q2,q3,q4,q5,q6,q7],q'));
    for n=1:cnt
       pos = Tg*[double(eePos(n,:)),1]';
       gripper(n,1:3) = pos(1:3,1)';
    end
    
end