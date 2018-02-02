function main(pitchStart,pitchInc,pitchEnd,rollStart,rollInc,rollEnd,yawStart,yawInc,yawEnd,xStart,xInc,xEnd,yStart,yInc,yEnd,zStart,zInc,zEnd)
%myFun - Description
%
% Syntax: main(dist,degStart,degEnd)
%
% Long description
f = 12.5;
dx = 11e-3;
dy = 11e-3;
W = 864;
H = 864;

M = [f/dx,0,W/2;0,f/dy,H/2;0,0,1];
%X = [-72.5, 17.5, 127.5, -72.5;18.75,18.75, 18.75,-56.25;1, 1, 1, 1];
%[x,y] = meshgrid(-90:90:90,-90:90:90);
%x = reshape(x,1,size(x,1)*size(x,2));
%y = reshape(y,1,size(y,1)*size(y,2));
x = [0 90 200 0 90];
y = [0 0 0 -75 -75];
X = [x;y];
X(3,:) = 1;

for i = xStart:xInc:xEnd
for j = yStart:yInc:yEnd
for k = zStart:zInc:zEnd
for l = pitchStart:pitchInc:pitchEnd
for m = rollStart:rollInc:rollEnd
for n = yawStart:yawInc:yawEnd
  Y = generateCoords(X,i,j,k,l,m,n,M);
  figure;
  scatter(Y(1,:),repmat(H,1,size(Y,2))-Y(2,:));
  xlim([0,W]);
  ylim([0,H]);
  
  disp(['x=',num2str(i),' y=',num2str(j),' z=',num2str(k),' pitch=',num2str(l),' roll=',num2str(m),' yaw=',num2str(n)]);
  [Angle,t] = mle(X,Y,M);
  disp(['x=',num2str(t(1)),' y=',num2str(t(2)),' z=',num2str(t(3)),' pitch=',num2str(Angle(1)),' roll=',num2str(Angle(2)),' yaw=',num2str(Angle(3))]);
end
end
end
end
end
end





% figure;
% plot(X,Y1,'LineWidth',1,'Marker','*');
% set(gca,'xlim',[-70,70]);
% set(gca,'ylim',[-70,70]);
% xlabel('\beta /��');
% ylabel('����\alpha /��');
% grid(gca,'on');
% grid(gca,'minor');
% legend('Z=40mm','Z=80mm','Z=120mm');
% print('C:\Users\Victo\OneDrive\��Ŀ\��Ŀ�Ӿ�\Y_alpha.png','-dpng','-r1200');

% figure;
% plot(X,Y2,'LineWidth',1,'Marker','*');
% set(gca,'xlim',[-70,70]);
% xlabel('\beta /��');
% ylabel('����\beta /��');
% grid(gca,'on');
% grid(gca,'minor');
% legend({'Z=40mm','Z=80mm','Z=120mm'},'Location','southeast');


% figure;
% plot(X,Y3,'LineWidth',1,'Marker','*');
% set(gca,'xlim',[-70,70]);
% xlabel('\beta /��');
% ylabel('����\gamma /��');
% grid(gca,'on');
% grid(gca,'minor');
% legend({'Z=40mm','Z=80mm','Z=120mm'},'Location','southeast');

% figure;
% plot(X,Y2-repmat(X,1,size(Y2,2)),'LineWidth',1,'Marker','*');
% set(gca,'xlim',[-70,70],'ylim',[-1.1*max(max(abs(Y2-repmat(X,1,size(Y2,2))))),1.1*max(max(abs(Y2-repmat(X,1,size(Y2,2)))))]);
% grid(gca,'on');
% grid(gca,'minor');
% xlabel('\beta /��');
% ylabel('E(\beta) /��');
% legend('Z=40mm','Z=80mm','Z=120mm');

end