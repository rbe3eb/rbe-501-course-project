function [] = jointSpaceMotionControl(qEqn,dqEqn,ddqEqn,HT,Jv,qi,dqi,ti,tf) 
    global eePos eePosErr torque tme;
    eePos = []; eePosErr = []; torque = []; tme = [];
    
    % Initial position and velocity error
    X0 = [qi; dqi];
    
    % Solve for joint position and velocity errors
    options = odeset('RelTol', 1e-4, 'AbsTol', [1e-4; 1e-4; 1e-4; 1e-4; 1e-4; 1e-4]);
    [T,X] = ode45(@(t,x)diffSolverJointSpace(t,x,qEqn,dqEqn,ddqEqn,HT,Jv,method),[ti tf],X0,options);
end

function dx = diffSolverJointSpace(t,x,qEqn,dqEqn,ddqEqn,HT,~)
	global eePos eePosErr torque tme;   
    syms tt;
    
    dx = zeros(14,1);
    Kp = 100*eye(7);
    Kv = 100*eye(7);
    
    % Compute current desired motion -> [position, velocity, acceleration]
    qDes = eval(subs(qEqn{1}, tt, t));
    dqDes = eval(subs(dqEqn{1}, tt, t));
    ddqDes = eval(subs(ddqEqn{1}, tt, t));
    
    % Obtain state variables -> [joint angle error, joint velocity error]
	curr_q = x(1:7,1);
	curr_dq = x(8:14,1);
    
    % Compute error in the position and velocity
    qErr = qDes - curr_q;
    dqErr = dqDes - curr_dq;
    
    % Compute end-effector position and end-effector position error
    curr_eePos = FK(HT{end}, curr_q);
    pDes = FK(HT{end}, qDes);
    curr_eePosErr = pDes - curr_eePos;
    
    % collect end-effector position and position error
    eePos = [eePos; curr_eePos'];
    eePosErr = [eePosErr; curr_eePosErr'];
    tme = [tme; t];
    
    % Collect dynamic model
    M = mMatrix(curr_q);
    C = cMatrix(curr_dq);
    G = gMatrix(curr_q);
    
    tau = M*(ddqDes + Kv*dqErr + Kp*qErr) + C.*curr_dq + G;

    % Update state vector
    dx(1:7,1) = curr_dq;
    dx(8:14,1) = inv(M)*(tau - C.*curr_dq - G);
    torque = [torque; tau'];
end