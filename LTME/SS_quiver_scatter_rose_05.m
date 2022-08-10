%%
st = 1;
en = length(depth);

% middle velocity
t = floor((depth-1)/2);
for i = 1:length(t)
    um(i) = u(t(i),i); 
    vm(i) = v(t(i),i);
%     up_m(i) = up(t(i),i);
end
% surface velocity
t = floor(depth)-4;
for i = 1:length(t)
    us(i) = u(t(i),i); 
    vs(i) = v(t(i),i);
%     up_s(i) = up(t(i),i);
end

%% 3개층 유속 quiver로 나타내기
figure()
subplot(311)
h=quiver(us, vs);
set(h,'showarrowhead','off');
ylabel('cm/s','fontsize',14)
text(100, -53,'Surface','fontsize',14);
ylim([-80 80]);
xlim([0 en]);
title('Sumjin ADCP bottom mooring data  3 layer velocity time series' );
set(gca,'xTick',[0:144:8928],'xTickLabel',[6:31,1:10],'xlim',[0 17*144]); 

subplot(312)
h=quiver(um, vm);
set(h,'showarrowhead','off');
ylabel('cm/s','fontsize',14)
text(100, -53,'Middle','fontsize',14);
ylim([-80 80]);
xlim([0 en])
set(gca,'xTick',[0:144:8928],'xTickLabel',[6:31,1:10],'xlim',[0 17*144]); 

subplot(313)
h=quiver(u(1,:), v(1,:));
set(h,'showarrowhead','off');
ylabel('cm/s','fontsize',14)
text(100, -53,'Bottom','fontsize',14);
ylim([-80 80]);
xlim([0 en]);
xlabel('time (day)')
set(gca,'xTick',[0:144:8928],'xTickLabel',[6:31,1:10],'xlim',[0 17*144]); 

%% 조류자료분석 - scatter plot, rose diagram
mu = nanmean(u,1);
mv = nanmean(v,1);
figure()
subplot(121)
plot(us,vs,'r.','MarkerSize',5)   % scatter plot 
hold on;
plot(um,vm,'g.','MarkerSize',5);hold on;
plot(u(1,:),v(1,:),'b.','MarkerSize',5);hold on;
plot([-100 100],[0 0],'k');plot([0 0],[-100 100],'k');% 중심표현
xlim([-100 100]); ylim([-100 100]);
legend('surface','middle','bottom');
axis equal
title('Scatter plot');xlabel('U comp (cm/sec)');ylabel('V comp (cm/sec)')

subplot(122)
polar(1,100);hold on;
plot(us,vs,'r.');
plot(um,vm,'g.');
plot(u(1,:),v(1,:),'b.');
title('Degrees Rose histogram')
degree = dir;


