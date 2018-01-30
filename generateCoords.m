function Y = generateCoords(X,x,y,z,pitch,roll,yaw,M)
%myFun - Description
%
% Syntax: Y = generateCoords(x,y,z,pitch,roll,yaw,M)
%
% Long description
X(3,:) = 0;
R = [cos(yaw/180*pi),0,sin(yaw/180*pi);0,1,0;-sin(yaw/180*pi),0,cos(yaw/180*pi)]*[1,0,0;0,cos(pitch/180*pi),-sin(pitch/180*pi);0,sin(pitch/180*pi),cos(pitch/180*pi)]*[cos(roll/180*pi),-sin(roll/180*pi),0;sin(roll/180*pi),cos(roll/180*pi),0;0,0,1];
Y = R*X + repmat([x;y;z],1,size(X,2));
Y(1,:) = Y(1,:)./Y(3,:);
Y(2,:) = Y(2,:)./Y(3,:);
Y(3,:) = 1;
Y = M*Y;

end