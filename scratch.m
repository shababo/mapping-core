figure;

offset = 100;

t = (0:size(traces,2)-1)/20000;

plot(t,mean(traces(33:38,:)),'r')
hold on;
plot(t,mean(traces(39:43,:))-offset,'g')
hold on;
plot(t,mean(traces(44:49,:)) - 2*offset,'c')
hold on;
plot(t,mean(traces(92:100,:)) - 3*offset,'b')
% hold on
% plot(t,mean(traces(142:146,:)) - 4*offset,'k')
hold off

legend({'950','900','850','800'})
title('Chronos Spectrum (one cell)')
xlabel('Time (sec)')
ylabel('Current (pA)')