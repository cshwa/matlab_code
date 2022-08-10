%% 유량과 조위 값
% 유량자료 가져오기 - 영산강홍수통제소 송정
% 1st column: 수위표수위(m)
% 2nd column: 유량(cms)
% 3rd column: 해발수위(m)
directory='D:\SilverStar\data\01_Sumjin_river';
filename='songjung201501.csv';
filename=fullfile(directory,filename);
riv01 = xlsread(filename);
riv01=riv01(:,2);

filename='songjung201502.csv';
filename=fullfile(directory,filename);
riv02 = xlsread(filename);
riv02=riv02(:,2);

filename='songjung201503.csv';
filename=fullfile(directory,filename);
riv03 = xlsread(filename);
riv03=riv03(:,2);

% 관측 기간 유량
a = riv02(2362:end);
b = riv03(1:782);
riv = [a' b']';

% 광양 조위값 : 1시간 간격 자료
% http://sms.khoa.go.kr/koofs/kor/observation/obs_past_search.asp?contents=1hour
% 시작: 2015 02 17 09 00
% 종료: 2015 03 06 10 00
tide = load('D:\SilverStar\data\01_Sumjin_tide\gyTG\gy_1h_2015_ing.txt');
tide = tide(1125:1547,6);

% 유량과 조위 값 그려보기
figure()
subplot(211)
l = length(riv);
x = [1:l];
plot(x,riv,'linewidth',2);
title('River discharge  (Feb. 17 - Mar. 6)','fontsize',14);
% xlabel('time (day)');
set(gca,'xTick',[0:144:8928],'xTickLabel',[17:28,1:7],'xlim',[0 17*144]);                
ylabel('cms','fontsize',14);

subplot(212)
plot(tide,'linewidth',2);
title('Tide (Feb. 17 - Mar. 6)','fontsize',14);
xlabel('Time (day)','fontsize',14);
ylabel('cm','fontsize',14) % left y-axis
set(gca,'xTick',[0:24:8928],'xTickLabel',[17:28,1:7],'xlim',[0 17*24]);                
