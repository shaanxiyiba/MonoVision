function M = crossMatrix(v)
%crossMatrix - Description
%
% Syntax: M = crossMatrix(v)
%
% Get cross product matrix from vector v
% Input:
% v - 3 x 1 vector
% M - 3 x 1 matrix

M = zeros(3,3);
M(1,2) = -v(3);
M(1,3) = v(2);
M(2,1) = v(3);
M(2,3) = -v(1);
M(3,1) = -v(2);
M(3,2) = v(1);

end