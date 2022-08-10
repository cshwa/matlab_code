
% 시간 간격: 10분
% 1st column: 수위(m)
% 2nd column: 유량
% 3rd column: 해발수위
%% river diacharge
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


%%
figure()
plot(riv(1:140),'k','linewidth',2);
hold on;
for i = 1: 4
    plot(NumDays(i), 0:0.1:150,'k-');
end
set(gca,'xtick',[1:10:150],'xticklabel',[1:10:150],'ylim',[0 150]);
ylabel('m^3/s','fontsize',12);
title('River Discharge','fontsize',14);


