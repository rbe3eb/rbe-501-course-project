% Compute joint values given end-effector position
function currQ = IK(T,Jv,Jw,pDes,rDes,currP,~,currQ)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    q = [q1; q2; q3; q4; q5; q6; q7];
    
    currQ = currQ';
    currR = rotMat(T,currQ);
    J = [Jv; Jw];
	n = 1;
    epsilon = 0.000001;
    err = 100;
    alpha = 1;
    while (n == 1 || norm((err)) > epsilon && n < 200)
        if isequal(mod(n,50),0)
            currQ = (pi).*rand(7,1);
        end
        % Evaluate jacobian given current joint values
        curr_J = eval(subs(J, q, currQ));
        
        % ERROR
        pErr = pDes - currP;
        rErr = rotm2Vector(rDes * currR');
        err = [pErr; rErr];
        
        % Update joint values
        currQ = currQ + alpha*pinv(curr_J)*(err);
        
        % Evaluate position given current joint values
        currT = eval(subs(T,[q1,q2,q3,q4,q5,q6,q7],currQ'));
        currR = currT(1:3,1:3);
        currP = currT(1:3,4);
        n = n + 1;
    end
end