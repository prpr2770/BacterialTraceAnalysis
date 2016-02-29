function getVelAcc_OfTracePair(traces)
% ---------------------------------------
% VEL-ACC plots

x1 = traces(:,1);
x2 = traces(:,2);

figure;
subplot(1,2,1)
plot(x1/norm(x1),'b')
hold on
v1 = diff(x1);
plot(v1/norm(v1),'g')
hold on
a1 = diff(v1);
plot(a1/norm(a1),'r')
hold off

subplot(2,1,2)
plot(x2/norm(x2),'b')
hold on
v2 = diff(x2);
plot(v2/norm(v2),'g')
hold on
a2 = diff(v2);
plot(a2/norm(a2),'r')
hold off
title('DISP, VEL and ACC')

end