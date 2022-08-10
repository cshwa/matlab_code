clc; clear all; close all;

%% data load - daily value : 12년부터 15년 현재 관측일 까지
pwd
directory = pwd;
filename = 'gy_dailyairtemp12_15.xlsx';
for i = 1:4
sheet = i;
xlRange = 'B2:M32';
% temp = xlsread(filename,sheet,xlRange);
[temp, txt, raw] = xlsread(filename,sheet,xlRange,'basic');

data = temp; 
[a b] = size(data);
data = reshape(data,[a*b,1]);
% at12 = snip(data,nan);
 eval(['at' num2str(i+11) '=snip(data,nan);']);
end
at13(at13 == 32768) = NaN;

for i = 5:8
sheet = i;
xlRange = 'B2:M32';
temp = xlsread(filename,sheet,xlRange,'basic');

data = temp; 
[a b] = size(data);
data = reshape(data,[a*b,1]);
% at12 = snip(data,nan);
 eval(['ys_at' num2str(i+7) '=snip(data,nan);']);
end


grey = [0.8,0.8,0.8];
FigHandle = figure;
set(FigHandle, 'Position', [100, 0, 600, 180]);
plot(at12,'.','color',grey);
hold on;

plot(at14,'.','color',grey);
plot(at13,'m.');
plot(at15,'k.-');
legend('12','13','14','15');

%--- 여수 기온 데이터 ----
% plot(ys_at12,'+','color',grey);
% hold on;
% plot(ys_at14,'+','color',grey);
% plot(ys_at13,'r+-');
% plot(ys_at15,'b+-');
%-------------------------
% 2015년 관측일
obs1 = juliandate(2015,02,26);
obs2 = juliandate(2015,03,06);
obs3 = juliandate(2015,05,13);
obs4 = juliandate(2015,05,19);
obs5 = juliandate(2015,08,16);
obs6 = juliandate(2015,08,27);

MoreDays = ['26-Feb-2015'; '03-Mar-2015'; '13-May-2015';'19-May-2015';'16-Aug-2015';'27-Aug-2015'];
NumDays = days365('01-jan-2015',MoreDays);

for i = 1:6
    plot(NumDays(i), -10:0.1:40,'b-','linewidth',1.5);
end
set(gca,'xtick',[1,32,60,91,121,151,182,212,243,273,304,334],'xticklabel',[1:12])
set(gca,'xlim',[1 365],'ylim',[-10 40],'fontsize',14);
ylabel('℃','fontsize',14);
title('Air Temperature','fontsize',14);

%% data load - daily precipitation : 12년부터 15년 현재 관측일 까지
pwd
directory = pwd;
filename = 'gy_precipitaion12_15.xlsx';
for i = 1:4
sheet = i;
xlRange = 'B2:M32';
temp = xlsread(filename,sheet,xlRange);
data = temp; 
[a b] = size(data);
data(isnan(data)) = 0;
data = reshape(data,[a*b,1]);
data(data == 32768) = NaN;
eval(['pr' num2str(i+11) '= snip(data,nan);']);
end
% at13(at13 == 32768) = NaN;

% 여수 데이터 처리용
% for i = 5:8
% sheet = i;
% xlRange = 'B2:M32';
% temp = xlsread(filename,sheet,xlRange);
% 
% data = temp; 
% [a b] = size(data);
% data = reshape(data,[a*b,1]);
% % at12 = snip(data,nan);
%  eval(['ys_at' num2str(i+7) '=snip(data,nan);']);
% end




%--- 여수 기온 데이터 ----
% plot(ys_at12,'+','color',grey);
% hold on;
% plot(ys_at14,'+','color',grey);
% plot(ys_at13,'r+-');
% plot(ys_at15,'b+-');
%-------------------------
% 2015년 관측일
obs1 = juliandate(2015,02,26);
obs2 = juliandate(2015,03,06);
obs3 = juliandate(2015,05,13);
obs4 = juliandate(2015,05,19);
obs5 = juliandate(2015,08,16);
obs6 = juliandate(2015,08,27);

MoreDays = ['26-Feb-2015'; '03-Mar-2015'; '13-May-2015';'19-May-2015';'16-Aug-2015';'27-Aug-2015'];
NumDays = days365('01-jan-2015',MoreDays);


grey = [0.8,0.8,0.8];
FigHandle = figure;
set(FigHandle, 'Position', [100, 0, 600, 180]);
plot(pr12,'.','color',grey);
hold on;
plot(pr13,'.','color',grey);
plot(pr14,'.','color',grey);
plot(pr15,'k.-');
legend('12','13','14','15');

for i = 1:6
    plot(NumDays(i), 00:0.1:200,'b-','linewidth',1.5);
end
set(gca,'xtick',[1,32,60,91,121,151,182,212,243,273,304,334],'xticklabel',[1:12])
set(gca,'xlim',[1 365],'ylim',[0 200],'fontsize',14);
ylabel('mm','fontsize',14);
title('Precipitation','fontsize',14);


%%
% filename = 'gy_aws.xlsx';
% filename = fullfile(directory,filename);
% data =  xlsread(filename);
% data = flipud(data);
% lowTemp = data(1:2:end,3);
% highTemp = data(1:2:end,4);
% lowHum = data(1:2:end,6);
% highHum = data(1:2:end,6);
% dayPreci = data(1:2:end,5);
% lowPress = data(1:2:end,8);
% highPress = data(1:2:end,9);

%{

filename = 'gy_dailymean_airtemp_past.xlsx';
filename = fullfile(directory,filename);
data =  xlsread(filename);
data = data(1:end-1,1:8);
[a b] = size(data);
airTemp = reshape(data, [a*b,1]);
airTemp = snip(airTemp,nan);

filename = 'gy_windspeed_ms_past.xlsx';
filename = fullfile(directory,filename);
data =  xlsread(filename);
data = data(1:end-1,1:8);
[a b] = size(data);
wind = reshape(data, [a*b,1]);
wind= snip(wind,nan);

filename = 'gy_precipitaion_mm_past.xlsx';
filename = fullfile(directory,filename);
data =  xlsread(filename);
data = data(1:end-1,1:8);
[a b] = size(data);
preci = reshape(data, [a*b,1]);
preci = preci([1:59,63:123,125:186,188:end]);

%}
%% 2015년 관측일
obs1 = juliandate(2015,02,26);
obs2 = juliandate(2015,03,06);
obs3 = juliandate(2015,05,13);
obs4 = juliandate(2015,05,19);
obs5 = juliandate(2015,08,16);
obs6 = juliandate(2015,08,27);

MoreDays = ['26-Feb-2015'; '03-Mar-2015'; '13-May-2015';'19-May-2015';'16-Aug-2015';'27-Aug-2015'];
NumDays = days365('01-jan-2015',MoreDays);
%% 그림 그리기
FigHandle = figure;
set(FigHandle, 'Position', [100, 0, 550, 400]);
subplot(311)
hold on;
plot(airTemp,'k.');
for i = 1:6
    plot(NumDays(i), -10:0.1:150,'b-','linewidth',1.5);
end
set(gca,'xtick',[1,32,60,91,121,151,182,212,243],'xticklabel',[1:9],'ylim',[-10 30],'fontsize',14);
ylabel('℃','fontsize',14);
title('Air Temperature','fontsize',20);

% subplot(313)
% plot(wind,'k.');
% hold on;
% for i = 1:6
%     plot(NumDays(i), -10:0.01:150,'b-');
% end
% set(gca,'xtick',[1,32,60,91,121,151,182,212,243],'xticklabel',[1:9],'ylim',[-10 30],'fontsize',14);
% set(gca,'ylim',[0 10],'fontsize',14);
% title('Wind speed','fontsize',24);
% ylabel('m/s','fontsize',14);

subplot(312)
plot(preci,'k.','linewidth',2);
hold on;
for i = 1:6
    plot(NumDays(i), -10:0.1:150,'b-','linewidth',1.5);
end
set(gca,'xtick',[1,32,60,91,121,151,182,212,243],'xticklabel',[1:9],'ylim',[-10 30],'fontsize',14);
set(gca, 'ylim',[0 100],'fontsize',14);
title('Daily Precipitation','fontsize',20);
ylabel('mm','fontsize',14);
% xlabel('time (month)','fontsize',18);



% save lowTemp.dat lowTemp -ascii
% save highTemp.dat highTemp -ascii
% save lowHum.dat lowHum -ascii
% save highHum.dat highHum -ascii

%% river diacharge
% 시간 간격: 10분
% 1st column: 수위(m)
% 2nd column: 유량
% 3rd column: 해발수위
directory = 'D:\SilverStar\data\01_Sumjin_river\songjung2015';
filename = 'songjung201501.csv';
filename = fullfile(directory,filename);
data =  xlsread(filename);
river = data;

filename = 'songjung201502.csv';
filename = fullfile(directory,filename);
data =  xlsread(filename);
river = [river; data];

filename = 'songjung201503.csv';
filename = fullfile(directory,filename);
data =  xlsread(filename);
river = [river; data];

filename = 'songjung201504.csv';
filename = fullfile(directory,filename);
data =  xlsread(filename);
river = [river; data];

filename = 'songjung201505.csv';
filename = fullfile(directory,filename);
data =  xlsread(filename);
river = [river; data];

for i = 1: length(river)/144
    riv(i,1) = mean(river(1+144*(i-1):144*i,2));
end



%% 2015년 관측일
obs1 = juliandate(2015,02,26);
obs2 = juliandate(2015,03,06);
obs3 = juliandate(2015,05,13);
obs4 = juliandate(2015,05,19);

MoreDays = ['26-feb-2015'; '03-mar-2015'; '13-may-2015';'19-may-2015'];
NumDays = days365('01-jan-2015',MoreDays);


%% 일평균값 그림 그리기
% figure()
% subplot(211)
% plot(airTemp(1:140),'k','linewidth',2);
% hold on;
% for i = 1: 4
%     plot(NumDays(i), -10:0.1:30,'k-');
% end
% set(gca,'xtick',[1:10:150],'xticklabel',[1:10:150],'ylim',[-10 30]);
% ylabel('℃','fontsize',12);
% text(3, 20,'Air Temperature','fontsize',14);

% subplot(414)
% plot((lowPress+highPress)/2,'k','linewidth',2);
% hold on;
% for i = 1: 4
%     plot(NumDays(i),1000:0.1:1038,'k-');
% end
% set(gca,'xtick',[1:10:150],'xticklabel',[1:10:150]);
% set(gca,'ylim',[1000 1038]);
% % text(3, 1005,'Air Pressure','fontsize',14);
% ylabel('hPa','fontsize',12);
% 
% subplot(413)
% plot(dayPreci,'k','linewidth',2);
% hold on;
% for i = 1: 4
%     plot(NumDays(i),1:100,'k-');
% end
% set(gca,'xtick',[1:10:150],'xticklabel',[1:10:150]);
% text(3, 70,'Daily Precipitation','fontsize',14);
% ylabel('mm','fontsize',12);



subplot(313)
plot(riv(1:140),'k','linewidth',2);
hold on;
for i = 1: 4
    plot(NumDays(i), 0:0.1:150,'k-');
end
set(gca,'xtick',[1:10:150],'xticklabel',[1:10:150],'ylim',[0 150]);
ylabel('m^3/s','fontsize',12);
% text(3, 100,'River Discharge','fontsize',14);
xlabel('Yearday','fontsize',12);