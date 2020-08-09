T7_ee = [1, 0, 0, 0
    0, -1, 0, 0
    0, 0, -1, -0.0615
    0, 0, 0, 1];

Tee_g = [1, 0, 0, 0
    0, 1, 0, 0
    0, 0, 1, 0.1250
    0, 0, 0, 1];

m = [0.500657, 0.9250];
centerOfMass_ee = {[-0.000281; -0.011402; -0.031080; 1], [0; 0; 0.0583; 1]};
centerOfMass_g = cell(1,2);

centerOfMass_g{1} = T7_ee*centerOfMass_ee{1};
centerOfMass_g{2} = Tee_g*centerOfMass_ee{2};

IT = inertiaTensor(centerOfMass_g,m);

function [IT] = inertiaTensor(centerOfMass, m)
    com = (centerOfMass{1} + centerOfMass{2}) / 2;
    x = com(1);
    y = com(2);
    z = com(3);
    M = m(1) + m(2);
    IT = struct;
    IT.xx = M*(y^2 + z^2);
    IT.yy = M*(x^2+z^2);
    IT.zz = M*(x^2+y^2);
    IT.xy = M*(x*y);
    IT.xz = M*(x*z);
    IT.zy = M*(y*z);
end