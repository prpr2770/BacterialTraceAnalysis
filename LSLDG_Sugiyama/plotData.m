function plotData(X,Y)

c=3;
[xx,yy]=meshgrid(-10:0.1:10);
gmmpdf=0.4*exp(-(yy).^2/2.0-(xx-c).^2/2.0)/(2*pi)...
    +0.3*exp(-(yy+c).^2/2.0-(xx+c).^2/2.0)/(2*pi)...
    +0.3*exp(-(yy-c).^2/2.0-(xx+c).^2/2.0)/(2*pi);

figure(1); subplot(121); contour(yy,xx,gmmpdf); hold on
figure(1); subplot(121); plot(X(1,:),X(2,:),'k.','MarkerSize',10); 
axis([-7 7 -7 7]); axis square; title('Input Data'); hold off;
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);

figure(1); subplot(122); contour(yy,xx,gmmpdf); hold on
figure(1); subplot(122); plot(Y(1,:),Y(2,:),'k.','MarkerSize',10); 
axis([-7 7 -7 7]); axis square; title('After Mode Seeking'); hold off;
set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
