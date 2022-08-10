% 수위자료 2000- 2005
% 데이터 이름: 2000.dat
% 첫번째 열 날짜, 1시간 간격 자료

% 자료 출처 사이트

clc; clear all; close all;

FigHandle = figure;
set(FigHandle, 'Position', [100, 0, 750, 700]);
xtic = cumsum(eomday(2013,1:12));
xtic_label = [1:12];


MoreDays = ['26-Feb-2015'; '03-Mar-2015'; '13-May-2015';'19-May-2015'; ...
    '16-Aug-2015';'27-Aug-2015';'21-Nov-2015';'28-Nov-2015'; ...
    '22-Mar-2015'; '14-Apr-2015'; ...
    '22-Jan-2015'; '03-Jun-2015'; '13-Aug-2015'; '30-Oct-2015';...
    '04-Nov-2015'; '15-Nov-2015'];
NumDays = days365('01-jan-2015',MoreDays);


for i = 2:5
%--- river ----------------------------------
eval(['load river' num2str(i+2011) '.txt']);
eval(['temp = river' num2str(i+2011) ]);

temp = temp(:,1:end);
[a b] = size(temp);
temp = reshape(temp,[a*b,1]);
temp(isnan(temp))=[]; % nan 제거하기
river(1:365,i-1) = temp(1:365,1);
%============================================
subplot(311)
hold on;
plot(river(1:365,i-1),'-','color',[0.15*i+0.12,0.15*i+0.12,0.15*i+0.12],'linewidth',1.5)
set(gca, 'xtick',[0 xtic],'xticklabel',xtic_label,'xlim',[0 364])

%--- temp -----------------------------------
eval(['load temp' num2str(i+2011) '.txt']);
eval(['data = temp' num2str(i+2011) ]);

data = data(:,1:end);
[a b] = size(data);
data = reshape(data,[a*b,1]);
data(isnan(data))=[]; % nan 제거하기
tempe(1:365,i-1) = data(1:365,1);
%============================================

subplot(312)
hold on;
plot(tempe(1:365,i-1),'-','color',[0.17*i,0.17*i,0.17*i],'linewidth',1.5)
set(gca, 'xtick',[0 xtic],'xticklabel',xtic_label,'xlim',[0 364])

%--- wind -----------------------------------
eval(['load wind' num2str(i+2011) '.txt']);
eval(['temp = wind' num2str(i+2011) ]);

temp = temp(:,1:end);
[a b] = size(temp);
temp = reshape(temp,[a*b,1]);
temp(isnan(temp))=[]; % nan 제거하기
wind(1:365,i-1) = temp(1:365,1);
%============================================
subplot(313)
hold on;
plot(wind(1:365,i-1),'-','color',[0.17*i,0.17*i,0.17*i],'linewidth',1.5)
set(gca, 'xtick',[0 xtic],'xticklabel',xtic_label,'xlim',[0 364])
end

subplot(311)
plot(mean(river,2),'b-','linewidth',2);
% legend('2012','2013','2014','2015','2016','Mean');
legend('2013','2014','2015','2016','Mean');
% for i = 1:16
%     plot(NumDays(i), -0:15:4000,'r.-','markersize',2);
% end
ylabel('[cms]','fontsize',14);
text(15, 3300,'River discharge','fontsize',20);


subplot(312)
plot(mean(tempe,2),'b-','linewidth',2);
% for i = 1:16
%     plot(NumDays(i), -10:0.1:40,'r.-','markersize',2);
% end
ylabel('[℃]','fontsize',14);
text(15, 33,'Temperature','fontsize',20);


subplot(313)
plot(mean(wind,2),'b-','linewidth',2);
% for i = 1:16
%     plot(NumDays(i), 0:0.05:10,'r.-','markersize',2);
% end
xlabel('Time (month)','fontsize',14);
text(15, 8.5,'Wind speed','fontsize',20);
ylabel('[m/s]','fontsize',14);
%%


% 2015년 관측일
% obs1 = juliandate(2015,02,26);
% obs2 = juliandate(2015,03,06);
% obs3 = juliandate(2015,05,13);
% obs4 = juliandate(2015,05,19);
% obs5 = juliandate(2015,08,16);
% obs6 = juliandate(2015,08,27);
% obs7 = juliandate(2015,11,21);
% obs8 = juliandate(2015,11,28);
% obs9 = juliandate(2015,03,22);
% obs10 = juliandate(2015,04,15);
% 
% MoreDays = ['26-Feb-2015'; '03-Mar-2015'; '13-May-2015';'19-May-2015'; ...
%     '16-Aug-2015';'27-Aug-2015';'21-Nov-2015';'28-Nov-2015'; ...
%     '22-Mar-2015'; '14-Apr-2015'];
% NumDays = days365('01-jan-2015',MoreDays);

% grey = [0.7,0.7,0.7];
% FigHandle = figure;
% set(FigHandle, 'Position', [100, 0, 600, 180]);
% plot(ss'.','color',grey);
% hold on;

subplot(311)
for i = 1:16
    plot(NumDays(i), 2000,'r*','markersize',10);
end
legend('2013','2014','2015','2016','Mean','CTD','location','north','Orientation','horizontal');
subplot(312)
for i = 1:16
    plot(NumDays(i), 20,'r*','markersize',10);
end

subplot(313)
for i = 1:16
    plot(NumDays(i), 6,'r*','markersize',10);
end


%}