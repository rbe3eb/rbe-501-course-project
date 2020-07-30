function dx = exampleHelperTimeBasedStateInputsPickPlace(t,x,motionModel,qEqn,dqEqn,ddqEqn,HT,Teg,Jv)
    syms tt q1 q2 q3 q4 q5 q6 q7 real;
    
    % Compute current desired motion -> [position, velocity, acceleration]
    qDes = eval(subs(qEqn, tt, t));
    dqDes = eval(subs(dqEqn, tt, t));
    ddqDes = eval(subs(ddqEqn, tt, t));
      
    targetState = [qDes; dqDes];
    
    curr_q = x(1:7,1);
	curr_dq = x(8:14,1);
    
    % Joint errors
    qErr = qDes - curr_q;
    dqErr = dqDes - curr_dq;
    
    state = [curr_q; curr_dq];
    
    % Compute state derivative
    dx = derivative(motionModel, state, targetState);
    disp(t);
end
