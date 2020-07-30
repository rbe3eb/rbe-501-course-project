% Compute joint values given end-effector position
function curr_q = IK(T,Jv,Jw,pDes,rDes,currP,currR,curr_q)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    q = [q1; q2; q3; q4; q5; q6; q7];
    
    curr_q = curr_q';
    %currR = eval(subs(T(1:3,1:3),[q1,q2,q3,q4,q5,q6,q7],curr_q'));
    J = [Jv; Jw];
	n = 1;
    epsilon = 0.00001;
    while ((norm([pDes; oDes] - [currP; currO])) > epsilon && n < 400)
        if isequal(mod(n,20), 0)
            curr_q = 0 + (pi-0).*rand(7,1);
            curr_p = FK(T,curr_q);
            currR = eval(subs(T(1:3,1:3),[q1,q2,q3,q4,q5,q6,q7],curr_q'));
            curr_o = rotm2Vector(currR);
            curr_op = [curr_p; curr_o];
        end
        % Evaluate jacobian given current joint values
        curr_J = eval(subs(J, q, curr_q));
        % Determine desired change in joint values
        pErr = currP - pDes;
        rErr = rDes * currR';
        
        qDelta = pinv(curr_J)*(opDes-curr_op);
        % Update joint values
        curr_q = curr_q + qDelta;
        % Evaluate position given current joint values
        curr_p = FK(T,curr_q);
        currR = eval(subs(T(1:3,1:3),[q1,q2,q3,q4,q5,q6,q7],curr_q'));
        curr_o = rotm2Vector(currR);
        curr_op = [curr_p; curr_o];
        n = n + 1;
    end
end