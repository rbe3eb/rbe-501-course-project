function [ phi ] = rotm2Vector(R)
    % Input: a rotation matrix C
    % Output: the rotational vector which describes the rotation C
    thetha = acos(0.5*(R(1,1) + R(2,2) + R(3,3) - 1));
    if (abs(thetha) < eps)
    	n = zeros(3,1);
    else
        n = 1/(2*sin(thetha))*[R(3,2) - R(2,3);
        R(1,3) - R(3,1);
        R(2,1) - R(1,2)];
    end
    phi = thetha*n;
end