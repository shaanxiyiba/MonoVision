function [R,t,err2] = orthogonalIteration(X,Y,R0,epsilon,maxIteration)
%orthogonalIteration - Description
%
% Syntax: [R,t] = orthognalIteration(X,Y,R0,epsilon)
%
% Orthogonal Iteration Algorithm for Pose Estimation
% Input:
% X - 3 x n world coordinates of n targets 
% Y - 3 x n focus unified image coordinates of n projections of targets 
% R0 - 3 x 3 initial rotation matrix
% epsilon - min error
% Output:
% R - 3 x 3 rotation matrix
% t - 3 x 1 translation vector
% err - indicate the transform error

n = size(X,2);
V = zeros(3,3,n);
G = zeros(3,9);
F = zeros(3,3);
X = X(1:3,:);
Q = zeros(3,n);
R = R0;

for i=1:n
  V(:,:,i) = Y(:,i)*Y(:,i)'/(Y(:,i)'*Y(:,i));
  F = F + (eye(3,3) - V(:,:,i));
  G = G + (eye(3,3) - V(:,:,i))*kron(X(:,i)',eye(3,3));
end
F = -inv(F);

i=0;
err = ones(2,5)*Inf;
while(true)

i = i+1;
t = F*G*reshape(R,9,1);

for j=1:n
  Q(:,j) = V(:,:,j)*(R*X(:,j)+t);
end

M = (Q-repmat(mean(Q,2),1,n))*(X - repmat(mean(X,2),1,n))';
[u,s,v] = svd(M);
R = v*diag([1,1,det(v*u')])*u';


pitch = -asin(R(2,3))/pi*180;
roll = atan(R(2,1)/R(2,2))/pi*180;
yaw = atan(R(1,3)/R(3,3))/pi*180;

err(mod(i,2)+1,mod((i-mod(i,2))/2,5)+1) = 0;
for j=1:n
  err(mod(i,2)+1,mod((i-mod(i,2))/2,5)+1) = err(mod(i,2)+1,mod((i-mod(i,2))/2,5)+1) + norm((eye(3,3)-V(:,:,j))*(R*X(:,j)+t),2);
end
err(mod(i,2)+1,mod((i-mod(i,2))/2,5)+1) = err(mod(i,2)+1,mod((i-mod(i,2))/2,5)+1)/n;

if((err(mod(i,2)+1,mod((i-mod(i,2))/2,5)+1) == min(min(err)))&& abs(max(err(mod(i,2)+1,:))-min(err(mod(i,2)+1,:)))<epsilon)
    err2 = min(min(err));
    break;
end

if(i>maxIteration)
    err2 = min(min(err));
    break;
end

end


end