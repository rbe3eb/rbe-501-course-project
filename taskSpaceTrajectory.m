function [pEqn,dpEqn,ddpEqn] = taskSpaceTrajectory(Xi,Xf,dXi,dXf,ti,tf)
    syms q1 q2 q3 q4 q5 q6 q7 real;
    
    % Trajectory polynomial
    pEqn = sym(zeros(3,1)); 
    dpEqn = sym(zeros(3,1)); 
    ddpEqn = sym(zeros(3,1));
    % Create polynomial equation
    for n=1:3
        [pEqn(n),dpEqn(n),ddpEqn(n)] = generatePoly(Xi(n),Xf(n),dXi(n),dXf(n),ti,tf);
    end
end

% Generate polynomial for trajectory
function [qEqn,dqEqn,ddqEqn] = generatePoly(qi,qf,dqi,dqf,ti,tf)
    syms tt a0 a1 a2 a3 real
    % Polynomial trajectory
    qEqn = a0 + (a1*tt) + (a2*tt^2) + (a3*tt^3);
    dqEqn = diff(qEqn, tt);
    ddqEqn = diff(dqEqn, tt);
    % Plug in initial conditions
    qEqn_i = subs(qEqn, tt, ti);
    dqEqn_i = subs(dqEqn, tt, ti);
    qEqn_f = subs(qEqn, tt, tf);
    dqEqn_f = subs(dqEqn, tt, tf);
    % Solve for coefficients
    eqn = [qi==qEqn_i, dqi==dqEqn_i, qf==qEqn_f, dqf==dqEqn_f];
    [c0, c1, c2, c3] = solve(eqn, [a0,a1,a2,a3]);
    % Final trajectory equations
    qEqn = subs(qEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
    dqEqn = subs(dqEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
    ddqEqn = subs(ddqEqn, [a0,a1,a2,a3], [c0,c1,c2,c3]);
end