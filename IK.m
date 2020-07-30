% Compute joint values given end-effector position
function curr_q = IK(T,Jv,Jw,pDes,oDes,curr_q)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    q = [q1; q2; q3; q4; q5; q6; q7];

    curr_q = curr_q';
    curr_p = FK(T,curr_q);
    J = [Jv; Jw];
    opDes = [pDes; oDes];
	n = 1;
    epsilon = 0.00001;
    while ((norm(opDes - curr_p)) > epsilon && n < 400)
        if isequal(mod(n,20), 0)
            curr_q = 0 + (pi-0).*rand(7,1);
            curr_p = FK(T,curr_q);
        end
        % Evaluate jacobian given current joint values
        curr_J = eval(subs(J, q, curr_q));
        % Determine desired change in joint values
        qDelta = pinv(curr_J)*(opDes-curr_p);
        % Update joint values
        curr_q = curr_q + qDelta;
        % Evaluate position given current joint values
        curr_p = FK(T,curr_q);
        n = n + 1;
    end
end