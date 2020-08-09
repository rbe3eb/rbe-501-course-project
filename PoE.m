function T = PoE()
    syms q1 q2 q3 q4 q5 q6 q7 real;
    HT = hTran();
    HT = simplify(HT{end},'Steps',50);
    
    qq = [q1; q2; q3; q4; q5; q6; q7];

    M = [1, 0, 0, 0
        0, 1, 0, -1.18/1000
        0, 0, 1, 1187.3/1000
        0, 0, 0, 1];
    
    v = cell(10,1);
    w = cell(10,1);
    q = cell(10,1);
    
    w{1} = [0; 0; -1];
    w{2} = [0; 1; 0];
    w{3} = [0; 0; -1];
    w{4} = [0; 1; 0];
    w{5} = [0; 0; -1];
    w{6} = [0; 1; 0];
    w{7} = [0; 0; -1];
    
    q{1} = [0; 0; 156.40]/1000;
    q{2} = [0; -5.4; 128.4]/1000;
    q{3} = [0; -6.4; 210.4]/1000;
    q{4} = [0; -6.4; 210.4]/1000;
    q{5} = [0; -6.4; 208.4]/1000;
    q{6} = [0; 0; 105.9]/1000;
    q{7} = [0; 0; 105.9]/1000;
    
    for n=1:7
       v{n} = cross(-w{n},q{n}); 
    end
    
    % Concatenated angular velocites
    Ws = [w{1}, w{2}, w{3}, w{4}, w{5}, w{6}, w{7}];
    % Concatenated linear velocities
    Vs = [v{1}, v{2}, v{3}, v{4}, v{5}, v{6}, v{7}];
    % Product of exponentials transformations
    ES = cell(7,1);
    % Derive the Forward Kinematics using Product of Exponentials
    % Approach
    I = eye(3);
    % build transformations
    T = eye(4);
    for i=1:7
        w = skew(Ws(:,i));
        th = qq(i);
        R = I+(sin(th)*w)+((1-cos(th))*(w^2));
        v = (I*th+(1-cos(th))*w+(th-sin(th))*w^2)*Vs(:,i);
        ES{i} = simplify([[R; 0, 0, 0], [v; 1]]);
        T = simplify(T * ES{i});
    end
    T = simplify(T*M,'Steps',50);  
end

function W = skew(w)
    W = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

function T = hTran()
	syms q1 q2 q3 q4 q5 q6 q7 real
    T = cell(7,1);
    
    % Base to frame 1
    T01 = [cos(q1), -sin(q1), 0, 0
           -sin(q1), -cos(q1), 0, 0
           0, 0, -1, 0.1564
           0, 0, 0, 1];
       
    % Frame 1 to frame 2
	T12 = [cos(q2), -sin(q2), 0, 0
        0, 0, -1, 0.0054
        sin(q2), cos(q2), 0, -0.1284
        0, 0, 0, 1];
    
    % Frame 2 to frame 3
	T23 = [cos(q3), -sin(q3), 0, 0
        0, 0, 1, -0.2104
        -sin(q3), -cos(q3), 0, -0.0064
        0, 0, 0, 1];
    
    % Frame 3 to frame 4
	T34 = [cos(q4), -sin(q4), 0, 0
        0, 0, -1, 0.0064
        sin(q4), cos(q4), 0, -0.2104
        0, 0, 0, 1];

    % Frame 4 to frame 5
	T45 = [cos(q5), -sin(q5), 0, 0
        0, 0, 1, -0.2084
        -sin(q5), -cos(q5), 0, -0.0064
        0, 0, 0, 1];
    
    % Frame 5 to frame 6
     T56 = [cos(q6), -sin(q6), 0, 0
         0, 0, -1, 0
         sin(q6), cos(q6), 0, -0.1059
         0, 0, 0, 1];
     
     % Frame 6 to frame 7
     T67 = [cos(q7), -sin(q7), 0, 0
         0, 0, 1, -0.1059
         -sin(q7), -cos(q7), 0, 0
         0, 0, 0, 1];
     
    T7_ee = [1, 0, 0, 0
        0, -1, 0, 0
        0, 0, -1, -0.0615
        0, 0, 0, 1];

    Tee_g = [1, 0, 0, 0
        0, 1, 0, 0
        0, 0, 1, 0.1250
        0, 0, 0, 1];

    % Frame 7 to gripper frame
    T6g = simplify(T67 * T7_ee * Tee_g);
    
    % Calculate homogeneous transformations
    T{1} = simplify(T01);
    T{2} = simplify(T01*T12);
    T{3} = simplify(T01*T12*T23);
    T{4} = simplify(T01*T12*T23*T34);
    T{5} = simplify(T01*T12*T23*T34*T45);
    T{6} = simplify(T01*T12*T23*T34*T45*T56);
    T{7} = simplify(T01*T12*T23*T34*T45*T56*T6g);
end