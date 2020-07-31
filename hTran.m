function [T,Teg,com] = hTran()
	syms q1 q2 q3 q4 q5 q6 q7 real
    
    T = cell(7,1); com = sym(zeros(3,7));
    % Link of center of mass wrt the precedent joint reference frame
	Lc = sym([[-0.000023; -0.010364; -0.073360; 1], ...
        [-0.000044; -0.099580; -0.013278; 1], ...
        [-0.000044; -0.006641; -0.117892; 1], ...
        [-0.000018; -0.075478; -0.015006; 1], ...
        [0.000001; -0.009432; -0.063883; 1], ...
        [0.000001; -0.045483; -0.009650; 1], ...
        [-0.0001; 0.0057; 0.0764; 1]]);
    
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
    
    T02 = simplify(expand(T01*T12));
    
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

    Teg = [1, 0, 0, 0
        0, 1, 0, 0
        0, 0, 1, 0.1250
        0, 0, 0, 1];
    
    % Calculate homogeneous transformations
    T{1} = simplify(T02);
    T{2} = simplify(T{1}*T23);
    T{3} = simplify(T{2}*T34);
    T{4} = simplify(T{3}*T45);
    T{5} = simplify(T{4}*T56);
    T{6} = simplify(T{5}*T67);
    T{7} = simplify(T{6}*T7_ee);
end
