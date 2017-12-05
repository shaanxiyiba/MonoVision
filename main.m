function main(dist,distStart,distEnd,degStart,degEnd,imageFormat)
%myFun - Description
%
% Syntax: main(dist,degStart,degEnd)
%
% Long description

for i = degStart:10:degEnd
for j = distStart:40:distEnd
disp(['dist=',num2str(j),' deg=',num2str(i)]);
I  = imread([getenv('MONO_DATA_PATH'),'\',num2str(j),'\',num2str(i),'.',imageFormat]);
X  = [degStart-360:10:degEnd-360]';
[alpha,T] = measure(I,0.2);
Y1((i-degStart)/10+1,(j-distStart)/40+1) = alpha(1,1);
Y2((i-degStart)/10+1,(j-distStart)/40+1) = alpha(2,1);
Y3((i-degStart)/10+1,(j-distStart)/40+1) = alpha(3,1);
end
end
figure;
plot(X,Y1,'LineWidth',1,'Marker','*');
set(gca,'xlim',[-70,70]);
set(gca,'ylim',[-70,70]);
xlabel('\beta /°');
ylabel('测量\alpha /°');
grid(gca,'on');
grid(gca,'minor');
legend('Z=40mm','Z=80mm','Z=120mm');
print('C:\Users\Victo\OneDrive\项目\单目视觉\Y_alpha.png','-dpng','-r1200');

figure;
plot(X,Y2,'LineWidth',1,'Marker','*');
set(gca,'xlim',[-70,70]);
xlabel('\beta /°');
ylabel('测量\beta /°');
grid(gca,'on');
grid(gca,'minor');
legend({'Z=40mm','Z=80mm','Z=120mm'},'Location','southeast');
print('C:\Users\Victo\OneDrive\项目\单目视觉\Y_beta.png','-dpng','-r1200');


figure;
plot(X,Y3,'LineWidth',1,'Marker','*');
set(gca,'xlim',[-70,70]);
xlabel('\beta /°');
ylabel('测量\gamma /°');
grid(gca,'on');
grid(gca,'minor');
legend({'Z=40mm','Z=80mm','Z=120mm'},'Location','southeast');
print('C:\Users\Victo\OneDrive\项目\单目视觉\Y_gamma.png','-dpng','-r1200');

figure;
plot(X,Y2-repmat(X,1,size(Y2,2)),'LineWidth',1,'Marker','*');
set(gca,'xlim',[-70,70],'ylim',[-1.1*max(max(abs(Y2-repmat(X,1,size(Y2,2))))),1.1*max(max(abs(Y2-repmat(X,1,size(Y2,2)))))]);
grid(gca,'on');
grid(gca,'minor');
xlabel('\beta /°');
ylabel('E(\beta) /°');
legend('Z=40mm','Z=80mm','Z=120mm');
print('C:\Users\Victo\OneDrive\项目\单目视觉\Y_error.png','-dpng','-r1200');

for i=1:3
RMSE(i) = sqrt(mean((Y2(:,i)-X).^2));    
disp(["RMSE ",num2str(i),"= ",num2str(RMSE(i))]);
end

end