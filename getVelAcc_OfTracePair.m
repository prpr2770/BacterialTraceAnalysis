function [vel, acc] = getVelAcc_OfTracePair(traces,Fs)
% ---------------------------------------
% Normalized VEL-ACC plots

x1 = traces(:,1);
x2 = traces(:,2);

figure;

subplot(3,2,1)
plot(x1/norm(x1),'b')
subplot(3,2,3)
v1 = Fs*diff(x1);
plot(v1/norm(v1),'g')
subplot(3,2,5)
a1 = Fs*diff(v1);
plot(a1/norm(a1),'r')

subplot(3,2,2)
plot(x2/norm(x2),'b')
subplot(3,2,4)
v2 = diff(x2);
plot(v2/norm(v2),'g')
subplot(3,2,6)
a2 = diff(v2);
plot(a2/norm(a2),'r')

set(gcf,'NextPlot','add');
axes;
h = title('Normalized DISP, VEL and ACC');
set(gca,'Visible','off');
set(h,'Visible','on');

vel = [v1 v2];
acc = [a1 a2];

end