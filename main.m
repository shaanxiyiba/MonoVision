function main(dist,distStart,distEnd,degStart,degEnd,imageFormat)
%myFun - Description
%
% Syntax: main(dist,degStart,degEnd)
%
% Long description

for i = degStart:10:degEnd
for j = distStart:40:distEnd
    
I  = imread([getenv('MONO_DATA_PATH'),'\',num2str(j),'\',num2str(i),'.',imageFormat]);
X  = [degStart-360:10:degEnd-360]';
Y((i-degStart)/10+1,(j-distStart)/40+1) = measure(I);
end
end
figure;
plot(X,Y,'LineWidth',1,'Marker','*');
set(gca,'xlim',[-100,100]);
set(gca,'ylim',[-100,100]);
xlabel('Desired rotation angle about Y axis\beta/бу');
ylabel('Measured rotation angle about Y axis\beta/бу');
grid(gca,'on');
grid(gca,'minor');

figure;
plot(X,Y-repmat(X,1,size(Y,2)),'LineWidth',1,'Marker','*');
set(gca,'xlim',[-90,90],'ylim',[-1.1*max(max(abs(Y-repmat(X,1,size(Y,2))))),1.1*max(max(abs(Y-repmat(X,1,size(Y,2)))))]);
grid(gca,'on');
grid(gca,'minor');
xlabel('Desired rotation angle about Y axis \beta/бу');
ylabel('Error of rotation angle measurement about Y axis \beta/бу');

end