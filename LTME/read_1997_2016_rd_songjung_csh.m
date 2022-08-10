% 수위자료 2000- 2005
% 데이터 이름: 2000.dat
% 첫번째 열 날짜, 1시간 간격 자료
% 화개
% clc; clear all; close all;

% csv read test
% 송정
clc; clear all; close all;



%% read multiple excel file


folder='D:\mepl\data\03_Sumjin_river\songjung\sj_00_16_year\';
filetype='*.xls';  % or xlsx
f=fullfile(folder,filetype);
d=dir(f);
for k=1:numel(d);
  data{k}=xlsread(fullfile(folder,d(k).name));
  temp = data{k};
  temp = temp(:,2:end);
  [a b] = size(temp);
  temp2 = reshape(temp, a*b,1);
  temp2 (isnan(temp2))=[]; % nan 열 제거
  to_rd(1:365,k) = temp2(1:365); % 1년을 365일 기준으로 함
  to_1d(1:365,k) = temp2(1:365)*24*3600; % 1년을 365일 기준으로 함
  eval(['RD' num2str(1996+k) '= temp2;'])   
end
%}
temp = cumsum(to_1d);
to_1y = mean(temp(end,:)); %1년 총량 20년 평균
temp = cumsum(to_1d(182:273,:));
to_3m = mean(temp(end,:));%7-9월 총량 평균
to_3m/to_1y


%% 
ss = reshape(to_rd,[],1);
mode(ss)
nanmean(ss)
nanmedian(ss)
% mean of dry season
d_rd = mean(mean(to_rd([1:180,273:365],:)));

% 7, 8, 9월의 양의 비율 구하기
t_rd = mean(sum(to_rd)); % 20년 평균 1년 총 유량
w_rd = mean(sum(to_rd(181:272,:))); % 7,8,9월 의 일 평균 유량
ra_rd = w_rd/t_rd;
%% 1997년부터 2016년 유량의 시간변화

figure()
for i = 1:20
    plot(to_rd(:,i),'color',[(1-i/40),i/40,i/100]);
    hold on;
end
tm_rd=sum(to_rd,2)/20;
plot(tm_rd,'k','linewidth',2);

xtic = cumsum(eomday(2013,1:12));
xtic_label = [1:12];
set(gca, 'xtick',[0 xtic],'xticklabel',xtic_label)
xlabel('Time (month)','fontsize',14);

%% River discharge
%{
in = 1; % interval
l = length(RD);
figure();
plot(1:l,RD(1:in:end)','b.');
legend('River discharge');
ylabel('cms','fontsize',14)

%% 유량 데이터 값 정리 
% 0 은 NAN 으로
RD(RD == 0 ) = NaN;
nanmean(RD)
nanmedian(RD)
mode(RD)


%%
in = 1; % interval
l = length(RD2);
figure();
plotyy(1:l,RD2(1:in:end),1:l,WL(1:in:end));
% plot(HL(1:in:end),'g');
legend('River discharge','Water level');
%%
% 2014년 하동과 송정의 수위, 유량 비교
directory='D:\SilverStar\data\01_Sumjin_river';
filename='songjung2014.xlsx';
filename=fullfile(directory,filename);
temp = xlsread(filename);
data = flipud(temp);
WL_s = data(:,1);
RD_s = data(:,2);
HL_s = data(:,3);

figure()
subplot(211)
plot(WL,'b');
hold on;
plot(WL_s,'k');
legend('Hadong','Songjung');
title('Water level (m)','fontsize',14);

subplot(212)
plot(RD2,'b');
hold on;
plot(RD_s,'k');
title('River discharge (cms)','fontsize',14);

%}