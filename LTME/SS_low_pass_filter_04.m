%% 36h low pass filter of depth mean velocity
% dum=load('pus9802.dat');
% t=dum(:,1); 
% x=dum(:,2);
x = nanmean(up);
del_t=1/6;
Nyq_fq=1/(2*del_t); cutoff_fq=1/36;
wc=cutoff_fq/Nyq_fq;
order=8;
[b,a]=butter(order,wc,'low');
x2=filtfilt(b,a,x);

figure()
hold on;
plot(x,'b');
plot(x2,'r');
plot(0:2453,0,'k-');
legend('depth mean raw data','36h low pass');
set(gca,'xTick',[0:144:8928],'xTickLabel',[6:31,1:10],'xlim',[0 17*144]); 
set(gca,'ylim',[-70 70]);
xlabel('Time (day)','fontsize',14);
ylabel('Current (cm/s)','fontsize',14);

