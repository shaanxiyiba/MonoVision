roll=0;
pitch =0;
yaw = 0;
t = [0,0,800]';
f = 12.5;      %% Focus length
dx = 5.5e-3;   %% Width of pixel
dy = 5.5e-3;   %% Height of pixel
M = [f/dx,0,640;0,f/dy,640;0,0,1];

R1 = [cos(yaw),0,sin(yaw);0,1,0;-sin(yaw),0,cos(yaw)];
R2 = [1,0,0;0,cos(pitch),-sin(pitch);0,sin(pitch),cos(pitch)];
R3 = [cos(roll),-sin(roll),0;sin(roll),cos(roll),0;0,0,1];

R = R1*R2*R3;

P = [-72.5, 17.5, 127.5, -72.5;18.75,18.75, 18.75,-56.25;0, 0, 0, 0];

p =R*P + repmat(t,1,4);

p = p./t(3);
p(3,:) = 1;

Y = M*p;

xx = Y(1,:);
yy = 1280 - Y(2,:);

figure;
scatter(xx,yy);



