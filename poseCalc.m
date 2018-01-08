clc;
clear all;
f = 12.5;
X = [0.0 90.0 200.0 0.0;0.0 0.0 0.0 -75.0];
Y = [452.029785 696.454285 958.67334 446.407532;633.687866 631.050354 628.339355 398.885437];
Y = (Y - repmat(640,2,4)).*5.5e-3;
K = zeros(8,8);
U = zeros(8,1);
C = zeros(8,1);
for i=1:4
  K(2*i-1:2*i,1:2) = X(1,i).*eye(2,2);
  K(2*i-1:2*i,3) = X(1,i).*[-Y(1,i)/f,-Y(2,i)/f]';
  K(2*i-1:2*i,4:5) = X(2,i).*eye(2,2);
  K(2*i-1:2*i,6:7) = eye(2,2);
  K(2*i-1:2*i,8) = [-Y(1,i)/f,-Y(2,i)/f]';
  C(2*i-1:2*i,1) = X(2,i).*[Y(1,i)/f,Y(2,i)/f]';
end

U = inv(K)*C;