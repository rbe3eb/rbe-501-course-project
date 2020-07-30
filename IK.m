% Compute joint values given end-effector position
function currQ = IK(T,Jv,Jw,pDes,rDes,currP,~,currQ)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    q = [q1; q2; q3; q4; q5; q6; q7];
    
    currQ = currQ';
    currR = rotMat(T,currQ);
    J = [Jv; Jw];
	n = 1;
    epsilon = 0.00001;
    err = 100;
    alpha = 0.5;
    while (n == 1 || norm((err)) > epsilon && n < 400)
        % Evaluate jacobian given current joint values
        curr_J = eval(subs(J, q, currQ));
        
        % ERROR
        pErr = pDes - currP;
        rErr = rotm2Vector(rDes * currR');
        err = [pErr; rErr];
        
        % Update joint values
        currQ = currQ + alpha*pinv(curr_J)*(err);
        
        % Evaluate position given current joint values
        currP = FK(T,currQ);
        currR = rotMat(T,currQ);
        n = n + 1;
    end
end