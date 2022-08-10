
% 여수 조위 관측소 2012-2015년도 값 그리기
grey = [0.8,0.8,0.8];
FigHandle = figure;
set(FigHandle, 'Position', [100, 0, 600, 180]);
hold on;

plot(ss(3,:),'.','color',grey);
plot(1:335,ss(1,1:335),'.','color',[0.87,0.87,0.87]);
plot(1:335,ss(2,1:335),'m.');
plot(1:231,ss(4,1:231),'k.');
legend('2012', '2014','2013','2015','Orientation','horizontal');
datetick('x', 'mm')
xlabel('Month','fontsize',14); 
ylabel('Temperature','fontsize',14);

MoreDays = ['26-Feb-2015'; '03-Mar-2015'; '13-May-2015';'19-May-2015';'16-Aug-2015';'27-Aug-2015'];
NumDays = days365('01-jan-2015',MoreDays);

for i = 1:6
    plot(NumDays(i), 0:0.1:30,'b-','linewidth',1.5);
end