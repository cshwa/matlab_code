% 수위자료 2000- 2005
% 데이터 이름: 2000.dat
% 첫번째 열 날짜, 1시간 간격 자료

% 자료 출처 사이트

clc; clear all; close all;
for i = 1:5
%--- river ----------------------------------
eval(['load river' num2str(i+2011) '.txt']);
eval(['temp = river' num2str(i+2011) ]);

temp = temp(:,1:end);
[a b] = size(temp);
temp = reshape(temp,[a*b,1]);
temp(isnan(temp))=[]; % nan 제거하기
river(1:365,i) = temp(1:365,1);
%============================================

%--- wind -----------------------------------
eval(['load wind' num2str(i+2011) '.txt']);
eval(['temp = wind' num2str(i+2011) ]);

temp = temp(:,1:end);
[a b] = size(temp);
temp = reshape(temp,[a*b,1]);
temp(isnan(temp))=[]; % nan 제거하기
wind(1:365,i) = temp(1:365,1);
%============================================

%--- temp -----------------------------------
eval(['load temp' num2str(i+2011) '.txt']);
eval(['data = temp' num2str(i+2011) ]);

data = data(:,1:end);
[a b] = size(data);
data = reshape(data,[a*b,1]);
data(isnan(data))=[]; % nan 제거하기
tempe(1:365,i) = data(1:365,1);
%============================================

end

%%
% 2015년 관측일
obs1 = juliandate(2015,02,26);
obs2 = juliandate(2015,03,06);
obs3 = juliandate(2015,05,13);
obs4 = juliandate(2015,05,19);
obs5 = juliandate(2015,08,16);
obs6 = juliandate(2015,08,27);
obs7 = juliandate(2015,11,21);
obs8 = juliandate(2015,11,28);

MoreDays = ['26-Feb-2015'; '03-Mar-2015'; '13-May-2015';'19-May-2015';'16-Aug-2015';'27-Aug-2015';'21-Nov-2015';'28-Nov-2015'];
NumDays = days365('01-jan-2015',MoreDays);

% grey = [0.7,0.7,0.7];
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 0, 600, 180]);
% plot(ss'.','color',grey);
% hold on;


%%
FigHandle = figure;
set(FigHandle, 'Position', [100, 0, 750, 550]);
subplot(311)
% plot(mean(river(:,1:3),2),'-','color',[0.7 0.7 0.7]);
area(mean(river(:,1:3),2),'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
% h.facecolor = [0.5 0.5 0.5];
hold on;
plot(river(:,4),'k.');
legend('Mean(2012-2014)','2015');
rm = max(mean(river(:,1:3),2));
for i = 1:8
    plot(NumDays(i), -0:15:rm,'b.-','markersize',2);
end
% set(gca,'xtick',[1:30:250*24*6],'xticklabel',[1:30:250*24*6]);
dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
dd = cumsum(dayspermonth);
set(gca,'xtick',dd(1:12),'xticklabel',[1:12],'ylim',[-10 30],'fontsize',14);
set(gca, 'ylim',[0 rm],'xlim',[0 365]);
title('Riverdischarge','fontsize',14);
ylabel('[cms]','fontsize',14);
% xlabel('Time (month)','fontsize',14);

subplot(313)
% plot(mean(wind(:,1:3),2),'-','color',[0.7 0.7 0.7]);
area(mean(wind(:,1:3),2),'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
hold on;
plot(wind(:,4),'k.');
wm = max(mean(wind(:,1:3),2));
legend('Mean(2012-2014)','2015');
for i = 1:8
    plot(NumDays(i), 0:0.1:8,'b.-','markersize',2);
end
% set(gca,'xtick',[1:30:250*24*6],'xticklabel',[1:30:250*24*6]);
dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
dd = cumsum(dayspermonth);
set(gca,'xtick',dd(1:12),'xticklabel',[1:12],'ylim',[-10 30],'fontsize',14);
set(gca, 'ylim',[0 8],'xlim',[0 365]);
title('Wind speed','fontsize',14);
ylabel('[m/s]','fontsize',14);

subplot(312)
% plot(mean(tempe(:,1:3),2),'-','color',[0.7 0.7 0.7]);
area(mean(tempe(:,1:3),2),'EdgeColor','none','FaceColor',[0.7 0.7 0.7]);
hold on;
plot(tempe(:,4),'k.');
legend('Mean(2012-2014)','2015');
wm = max(mean(tempe(:,1:3),2))
for i = 1:8
    plot(NumDays(i), -5:0.5:35,'b.-','markersize',2);
end
% set(gca,'xtick',[1:30:250*24*6],'xticklabel',[1:30:250*24*6]);
dayspermonth = [0 31 28 31 30 31 30 31 31 30 31 30 31];
dd = cumsum(dayspermonth);
set(gca,'xtick',dd(1:12),'xticklabel',[1:12],'ylim',[-10 30],'fontsize',14);
set(gca, 'ylim',[-5 35],'xlim',[0 365]);
title('Temperature','fontsize',14);
ylabel('[℃]','fontsize',14);
xlabel('Time (month)','fontsize',14);

