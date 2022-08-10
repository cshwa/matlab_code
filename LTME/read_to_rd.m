load to_rd.mat
%%
%--- 원하는 년도: sec_yr ---------------
sec_yr =[1:17];
months = {  'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 800, 300]);
hold on;
plot(1:365,mean(to_rd,2),'k-','linewidth',1.5);
plot(1:365,to_rd(:,sec_yr),'-','linewidth',1,'color',[0.8 0.8 0.8]);
plot(1:365,mean(to_rd,2),'k-','linewidth',1.5);
set(gca,'xtick',[0 cumsum(eomday(2001,1:11))],'xticklabel',months,'fontsize',12);
ylabel('유량 (m^3s^-^1)','fontsize',14);
lh=legend('18-yr Mean','2000-2017');
set(lh,'fontsize',14);
xlim([0 365]);
%--- save figure --------------------------------------
outpath = '.\figure\';
out=[outpath 'discharge_2000_2017'];
set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=8;
y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,out,'tif');
%------------------------------------------------------

%%
%--- 2004년 관측일 표시 ----------------------------------------------------
% 대조: 2004.10.16
ch_data = to_rd(:,5);
plot(1:365, ch_data,'b','linewidth',1);
plot(290, ch_data(290),'bo','MarkerFaceColor','b','markersize',5);
%--- 2005년 관측일 표시 ----------------------------------------------------
% 대조: 2005.01.29,소조:2005.01.19,2005.04.16
ch_data = to_rd(:,6);
plot(1:365, ch_data,'r','linewidth',1);
plot(29, ch_data(29),'ro','MarkerFaceColor','r','markersize',5);
plot(19, ch_data(19),'ro','MarkerFaceColor','r','markersize',5);
plot(106, ch_data(106),'ro','MarkerFaceColor','r','markersize',5);
%--- 2006년 관측일 표시 ----------------------------------------------------
% 대조: 2006.03.30, 08.10, 11.06, 소조: 2006.01.24,08.03, 10.31
ch_data = to_rd(:,7);
plot(1:365, ch_data,'k','linewidth',1);
plot(24, ch_data(24),'ko','MarkerFaceColor','k','markersize',5);
plot(89, ch_data(89),'ko','MarkerFaceColor','k','markersize',5);
plot(222, ch_data(222),'ko','MarkerFaceColor','k','markersize',5);
plot(215, ch_data(215),'ko','MarkerFaceColor','k','markersize',5);
plot(304, ch_data(304),'ko','MarkerFaceColor','k','markersize',5);
plot(310, ch_data(310),'ko','MarkerFaceColor','k','markersize',5);
%--- 2013년 관측일 표시 ----------------------------------------------------
ch_data = to_rd(:,17);
plot(1:365, ch_data,'g','linewidth',1);
plot(303, ch_data(303),'go','MarkerFaceColor','g','markersize',5);
plot(309, ch_data(309),'go','MarkerFaceColor','g','markersize',5);
set(gca,'xtick',[0 cumsum(eomday(2001,1:11))],'xticklabel',months);
axis([0 365 0 4300])

%--- save figure --------------------------------------
outpath = '.\figure\';
out=[outpath,'River_','Observation'];
set(gcf,'renderer','painter');
set(gcf, 'PaperUnits', 'inches');
x_width=12;
y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
saveas(gcf,out,'tif');
%------------------------------------------------------
%%
% 20년 평균
mean(mean(to_rd))

% 드라이 시즌 평균 : 11,12,1,2 월
d1=sum(eomday(2001,1:10)); % 1월부터 10월까지 일수
d2=sum(eomday(2001,1:2)); % 2월까지 일 수 

mean(mean(to_rd([1:d2,d1+1:end],:)))
% 웻 시즌 평균: 7,8,9,10
w1=sum(eomday(2001,1:6)); % 1월부터 10월까지 일수
w2=sum(eomday(2001,1:9)); % 2월까지 일 수 

mean(mean(to_rd(w1+1:w2,:)))