% 섬진강 adcp mooring
% N34 59 13.3 
% E127 46 31.8
% 관측간격 1m
% First cell 1.5m
% 관측 시간 간격: 10 min
% 사용한 자료: 2015-2-17-06:00 - 2015-3-06-16:19

clc; clear all; close all;

% mat 형식 데이터 load
load all.mat

st = 1; en = 2511;
pg1 = SerPG4([st:en],:);

% pg1으로 관측 시작/끝 시점 찾기
% 100과 2400은 시작과 끝 지점에서 어느정도 떨어진 값으로 임의로 잡음
st = min(find(pg1(1:100,1) >= 99));
en = min(find(pg1(2400:end,1) < 99))+2400-2;

depth = AnDepthmm([st:en],:); % unit: mm
mag = SerMagmmpersec([st:en],:);
dir = SerDir10thDeg([st:en],:)/10;
v = SerNmmpersec([st:en],:);
u = SerEmmpersec([st:en],:);
w = SerVmmpersec([st:en],:);
pg1 = SerPG4([st:en],:);
pg4 = SerPG4([st:en],:);
error = SerErmmpersec([st:en],:);
coabg = SerCAcnt([st:en],:);
echoavg = SerEAAcnt([st:en],:);
dd = depth/1000;
depth = floor(depth/1000);  % unit: m

save dd.txt dd -ascii


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
title('Songjung - 2015 winter (Feb. to Mar.)');
% xlabel('time (day)');
ylabel('river discharge (cms)') % left y-axis
set(gca,'xTick',[0:144:8928],'xTickLabel',[0:1:31],'xlim',[0 17*144]);                
ylabel('cms');

subplot(212)
plot(tide,'linewidth',2);
title('Songjung - 2015 winter (Feb. to Mar.)');
xlabel('time (day)');
ylabel('tide') % left y-axis
set(gca,'xTick',[0:24:8928],'xTickLabel',[0:1:17],'xlim',[0 17*24]);                
ylabel('cm');



%% test : u,v 분해와 검증
% uu = mag.*sind(dir);
% vv = mag.*cosd(dir);
% speed = sqrt(uu.^2+vv.^2);
% direction = -atand(uu./vv);
% d2 = atand(u./v);
% k = dir-direction;
u_org = u;
u = u_org;
%% depth 자료로 cutting
for i = 1:length(depth)
    u(i,depth(i)-1:end) = NaN;
    v(i,depth(i)-1:end) = NaN;
    mag(i,depth(i)-1:end) = NaN;      
end

%% shadowzone 판단용
% for i = 1:length(depth)
%     u(i,depth(i)+4:20) = NaN;
%     v(i,depth(i)+4:20) = NaN;
%     mag(i,depth(i)+4:20) = NaN;  
% end

u = flipud(rot90(u(:,1:15)/10));
v = flipud(rot90(v(:,1:15)/10));
mag = flipud(rot90(mag(:,1:end)/10));

mag(mag == -32768) = nan;
mean_u = nanmean(u);
mean_v = nanmean(v);


%% rotate angle
angleV = -45;
theta = dir;
%Magnitude remains constant.  Convert to new theta

theta2 = theta - angleV;
%Compensate if new angle has passed 0 degrees
over_ang = find(theta2 > 360);
if isempty(over_ang) == 0;
    theta2(over_ang) = theta2(over_ang)-360;
end

%Reconfigure to new coordinate system
up = mag.*sind(theta2');
across = mag.*cosd(theta2');

% figure()
% mean_up = nanmean(up);
% plot(mean_up);
% hold on;
% plot([0:30:4531],0,'k.','markersize',2)
% set(gca,'xTick',[0:144:4752]);                
% set(gca,'xTickLabel',[0:1:33]);
% xlim([0 4531]); %2014년 관측과 기간을 맞춰봄
% ylim([-150 100]);

figure('name','check:ratate angle');
subplot(211);
contourf(across); shading flat;
ylabel('across');
% hold on; plot(dd,'linewidth',2)
caxis([-50 50]);colorbar
set(gca,'xTick',[0:144:8928],'xTickLabel',[0:1:31],'xlim',[0 17*144]); 
set(gca,'yTick',[-2:2:19.5],'yTickLabel',[0:2:20],'ylim',[-1 14]); 

subplot(212);
contourf(up); shading flat;
ylabel('up');
caxis([-80 80]);colorbar
set(gca,'xTick',[0:144:8928],'xTickLabel',[0:1:31],'xlim',[0 17*144]); 
set(gca,'yTick',[-1:2:15],'yTickLabel',[0:2:20],'ylim',[-1 14]); 
% hold on; plot(depth,'linewidth',2)


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
set(gca,'xTick',[0:144:8928],'xTickLabel',[0:1:31],'xlim',[0 17*144]); 
xlabel('time (day)');
%%
st = 1;
en = length(depth);

% figure()
% 
% subplot(511)
% x = [0:4533];
% [hAx,hLine1,hLine2] = plotyy(x,riv,x,level);
% title('Songjung river discharge - August 2014','fontsize',14);
% xlabel('time (day)');
% ylabel(hAx(1),'river discharge (cms)') % left y-axis
% ylabel(hAx(2),'water level (m)') % right y-axis
% set(hAx(1),'xTick',[0:144:8928],'xTickLabel',[1:1:31],'xlim',[0 en]);                
% set(hAx(2),'xTick',[0:144:8928],'xTickLabel',[1:1:31],'xlim',[0 en]);
% ylabel('cms');
% colorbar
% 
% subplot(512)
% plot(tide);
% en = length(tide);
% % ylim([-2 18]);
% ylabel('tide (cm)');
% % set(gca,'yTick',[-0.5:5:19.5],'yTickLabel',[0:5:20],'ylim',[0 15]); 
% set(gca,'xTick',[0:24:4700],'xlim',[0 en]);                
% set(gca,'xTickLabel',[1:1:31])
% % caxis([-60,60])
% colorbar
% % hold on;
% % plot(dd,'linewidth',2)
% en = length(depth);
% 
% subplot(513)
% contourf(mag(1:16, st:en)/10);shading flat;
% title('Sumjin river ADCP mooring - 2014-8-01-00:00 - 2014-09-16:00','fontsize',14);
% % title('Sumjin river ADCP mooring - 2014-7-30-06:00 - 2014-09-16:00');
% ylim([-2 18]);ylabel('magnitude');
% set(gca,'yTick',[-0.5:5:19.5],'yTickLabel',[0:5:20],'ylim',[0 15]); 
% set(gca,'xTick',[0:144:4700],'xlim',[0 en]);                
% set(gca,'xTickLabel',[1:1:31])
% % caxis([-60,60])
% colorbar
% % hold on;
% % plot(dd,'linewidth',2)
% 
% subplot(514)
% contourf(u(1:16, st:en)/10);shading flat;
% ylim([-2 18]);
% ylabel('u');
% set(gca,'yTick',[-0.5:5:19.5],'yTickLabel',[0:5:20],'ylim',[0 15]);                
% set(gca,'xTick',[0:144:4700],'xlim',[0 en]);                
% set(gca,'xTickLabel',[1:1:31])
% caxis([-100,100])
% colorbar
% 
% subplot(515)
% contourf(v(1:16, st:en)/10);shading flat;
% ylim([-2 18]); ylabel('v');
% set(gca,'yTick',[-0.5:5:19.5],'yTickLabel',[0:5:20],'ylim',[0 15]); 
% set(gca,'xTick',[0:144:4700],'xlim',[0 en]);                
% set(gca,'xTickLabel',[1:1:31]);
% caxis([-100,100])
% colorbar

%%
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


figure()
subplot(311)
h=quiver(us, vs);
set(h,'showarrowhead','off');
ylabel('surface (cm/s)')
ylim([-80 80]);xlim([0 4320]);
title('Sumjin ADCP bottom mooring data  3 layer velocity time series' );
set(gca,'xTick',[0:144:4700]);                
set(gca,'xTickLabel',[0:1:33])
subplot(312)
h=quiver(um, vm);
set(h,'showarrowhead','off');
ylabel('middle (cm/s)')
ylim([-80 80]);xlim([0 4320])
set(gca,'xTick',[0:144:4700]);                
set(gca,'xTickLabel',[0:1:33])
subplot(313)
h=quiver(u(1,:), v(1,:));
set(h,'showarrowhead','off');
ylabel('bottom (cm/s)')
ylim([-80 80]);
xlim([0 4320]);
xlabel('time (day)')
set(gca,'xTick',[0:144:4700]);                
set(gca,'xTickLabel',[0:1:33]);
%% 조류자료분석
mu = nanmean(u,1);
mv = nanmean(v,1);
figure()
subplot(121)
plot(us,vs,'r.','MarkerSize',5)   % scatter plot 
hold on;
plot(um,vm,'g.','MarkerSize',5);hold on;
plot(u(1,:),v(1,:),'b.','MarkerSize',5);hold on;
plot([-100 100],[0 0],'k');plot([0 0],[-100 100],'k');% 중심표현
xlim([-200 200]); ylim([-200 200]);
legend('surface','middle','bottom');
axis equal
title('Scatter plot');xlabel('U comp (cm/sec)');ylabel('V comp (cm/sec)')

subplot(122)
polar(1,100);hold on;
plot(mu,mv,'b.');
title('Degrees Rose histogram')

degree = dir;

%% 정북기준 유향의 histogram
% figure()
% hist(degree,36);   % 30도 간격으로 hisgogram 즉 360을 12로 나눈 각도
% set(gca,'xTick',[0:10:360]);      
%% 3개층의 u,v 유속크기
% figure()
% usl = length(us);
% x = [1: usl];
% subplot(311)
% [hAx,hLine1,hLine2] = plotyy(x,us,x,vs);
% text(2,80,'surface');
% ylabel(hAx(1),'u(cm/s)') % left y-axis
% ylabel(hAx(2),'v(cm/s)') % right y-axis
% set(hAx(1),'xTick',[0:144:4700],'xTickLabel',[0:1:33],'ylim',[-130 100],'xlim',[0 4530]);                
% set(hAx(2),'xTick',[0:144:4700],'xTickLabel',[0:1:33],'ylim',[-130 100],'xlim',[0 4530]); 
% 
% % set(gca,'xTick',[0:18:4700],'xTickLabel',[0:3:33]);  
% subplot(312)
% [hAx,hLine1,hLine2] = plotyy(x,um,x,vm);
% text(2,80,'middle');
% ylabel(hAx(1),'u(cm/s)') % left y-axis
% ylabel(hAx(2),'v(cm/s)') % right y-axis
% set(hAx(1),'xTick',[0:144:4700],'xTickLabel',[0:1:33],'ylim',[-130 100],'xlim',[0 4530]);                
% set(hAx(2),'xTick',[0:144:4700],'xTickLabel',[0:1:33],'ylim',[-130 100],'xlim',[0 4530]); 
% 
% subplot(313)
% ub = u(1,:); 
% vb = v(1,:);
% [hAx,hLine1,hLine2] =plotyy(x,ub,x,vb);
% text(2,80,'bottom');
% ylabel(hAx(1),'u(cm/s)') % left y-axis
% ylabel(hAx(2),'v(cm/s)') % right y-axis
% set(hAx(1),'xTick',[0:144:4700],'xTickLabel',[0:1:33],'ylim',[-130 100],'xlim',[0 4530]);                
% set(hAx(2),'xTick',[0:144:4700],'xTickLabel',[0:1:33],'ylim',[-130 100],'xlim',[0 4530]); 


%% 일주일평균치를 구하기 위한 변수
% u1 = u(1:13,145:144*8);
% v1 = v(1:13,145:144*8);
% u2 = u(1:13,144*8+1:144*15);
% v2 = v(1:13,144*8+1:144*15);
% u3 = u(1:13,144*15+1:144*22);
% v3 = v(1:13,144*15+1:144*22);
% u4 = u(1:13,144*22+1:144*29);
% v4 = v(1:13,144*22+1:144*29);

% neap1 = up(1:13,145:144*8);
% spring1 = up(1:13,144*8+1:144*15);
% neap2 = up(1:13,144*15+1:144*22);
% spring2 = up(1:13,144*22+1:144*29);
% 
% l = length(neap1);
% for i=1:1
% neap1(:,i) = neap1(:,i) - nanmean(neap1(:,i));
% spring1(:,i) = spring1(:,i) - nanmean(spring1(:,i));
% neap2(:,i) = neap2(:,i) - nanmean(neap2(:,i));
% spring2(:,i) = spring2(:,i) - nanmean(spring2(:,i));
% end
% 
% figure()
% subplot(411)
% contourf(neap1); shading flat;
% ylabel('neap1');
% colorbar;
% subplot(412)
% contourf(spring1); shading flat;
% ylabel('spring1');
% colorbar;
% subplot(413)
% contourf(neap2); shading flat;
% ylabel('neap1');
% colorbar;
% subplot(414)
% contourf(spring2); shading flat;
% ylabel('spring1');
% colorbar;

%% progressive vector 표현 - 전체기간
% 3개층에 대한 조류 연속벡터도를 작성
% input data: TC_data.dat
% m-file: progvec_exam2.m

% clc;clear all;close all;
% Data 불러오기
% [year,mm,dd,hr,min,u_up,v_up,u_mid,v_mid,...
%    u_down,v_down]=textread('TC_data.dat','');
ti = 600;          % time interval 1hour = 3600sec
n = length(dd);     % data의 size

% 시작점 위치 (0,0)
% 표층 data set 작성
nu_up = zeros(n,1); nv_up =zeros(n,1);
for i=2:1:n
    nu_up(i) = nu_up(i-1) + us(i-1)/100000*ti;
    nv_up(i) = nv_up(i-1) + vs(i-1)/100000*ti;
end
% 중층 data set 작성
nu_mid = zeros(n,1); nv_mid =zeros(n,1);
for i=2:1:n
    nu_mid(i) = nu_mid(i-1) + um(i-1)/100000*ti; 
    nv_mid(i) = nv_mid(i-1) + vm(i-1)/100000*ti;  
end
% 저층 data set 작성
nu_down = zeros(n,1); nv_down =zeros(n,1);
for i=2:1:n
    nu_down(i) = nu_down(i-1) + u(1,i-1)/100000*ti;
    nv_down(i) = nv_down(i-1) + v(1,i-1)/100000*ti;
end

figure()
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nu_up,nv_up,'r');hold on;     % pgrogressive vector 도시
plot(nu_mid,nv_mid,'g');hold on;
plot(nu_down,nv_down,'b');hold on;

% 기준점에서 연결한 점선 
plot([-100 100],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-100 120],':k');
axis([-100,100 -100,120]);
% axis equal
% axis auto
legend('start','Upper layer','middle layer','bottom layer',...
      'Location','Southeast');
title('Progressive Vector Diagram at 3 Layer');          % 제목
xlabel('x-direction(km)'); ylabel('y-direction(km)');     % y축 이름

%% progressive vector 표현 - 전체기간
% 3개층에 대한 조류 연속벡터도를 작성
% input data: TC_data.dat
% m-file: progvec_exam2.m

% clc;clear all;close all;
% Data 불러오기
% [year,mm,dd,hr,min,u_up,v_up,u_mid,v_mid,...
%    u_down,v_down]=textread('TC_data.dat','');
ti = 600;          % time interval 1hour = 3600sec
n = length(dd);     % data의 size
ub = u(1,:); vb = v(1,:);
stn = 145; % neap 시기 유량 적을때 16일~18일
enn = 1152; 
sts = 1153; % spring 시기 유량 적을때
ens = 2160;
us_n = us(stn:enn); vs_n = vs(stn:enn);
um_n = um(stn:enn); vm_n = vm(stn:enn);
ub_n = ub(stn:enn); vb_n = vb(stn:enn);
us_s = us(sts:ens); vs_s = vs(sts:ens);
um_s = um(sts:ens); vm_s = vm(sts:ens);
ub_s = ub(sts:ens); vb_s = vb(sts:ens);


% 시작점 위치 (0,0)
% 표층 data set 작성
n = length(us_n);
nun_up = zeros(n,1); nvn_up =zeros(n,1);
for i=2:1:n
    nun_up(i) = nun_up(i-1) + us_n(i-1)/100000*ti;
    nvn_up(i) = nvn_up(i-1) + vs_n(i-1)/100000*ti;
end

% 중층 data set 작성
nun_mid = zeros(n,1); nvn_mid =zeros(n,1);
for i=2:1:n
    nun_mid(i) = nun_mid(i-1) + um_n(i-1)/100000*ti; 
    nvn_mid(i) = nvn_mid(i-1) + vm_n(i-1)/100000*ti;  
end
% 저층 data set 작성
nun_down = zeros(n,1); nvn_down =zeros(n,1);
for i=2:1:n
    nun_down(i) = nun_down(i-1) + ub_n(i-1)/100000*ti;
    nvn_down(i) = nvn_down(i-1) + vb_n(i-1)/100000*ti;
end
% 시작점 위치 (0,0)
% 표층 data set 작성
nus_up = zeros(n,1); nvs_up =zeros(n,1);
for i=2:1:n
    nus_up(i) = nus_up(i-1) + us_s(i-1)/100000*ti;
    nvs_up(i) = nvs_up(i-1) + vs_s(i-1)/100000*ti;
end

% 중층 data set 작성
nus_mid = zeros(n,1); nvs_mid =zeros(n,1);
for i=2:1:n
    nus_mid(i) = nus_mid(i-1) + um_s(i-1)/100000*ti; 
    nvs_mid(i) = nvs_mid(i-1) + vm_s(i-1)/100000*ti;  
end
% 저층 data set 작성
nus_down = zeros(n,1); nvs_down =zeros(n,1);
for i=2:1:n
    nus_down(i) = nus_down(i-1) + ub_s(i-1)/100000*ti;
    nvs_down(i) = nvs_down(i-1) + vb_s(i-1)/100000*ti;
end

figure('name','weekly variation')
set(gcf,'units', 'pixels', 'pos',[100 100 850 350])
subplot(141)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nun_up,nvn_up,'r');hold on;     % pgrogressive vector 도시
plot(nun_mid,nvn_mid,'g');hold on;
plot(nun_down,nvn_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
xlim([-100,100]);
ylim([-200,100]);
text(-70, 50,'Neap','fontsize',14);
axis equal
% legend('start','Upper layer','middle layer','bottom layer',...
%       'Location','northEast');
title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');ylabel('y-direction(km)');     % y축 이름

subplot(142)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nus_up,nvs_up,'r');hold on;     % pgrogressive vector 도시
plot(nus_mid,nvs_mid,'g');hold on;
plot(nus_down,nvs_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
text(-70, 50,'Spring','fontsize',14);
xlim([-100,100]);
ylim([-200,100]);
axis equal
legend('start','Upper layer','middle layer','bottom layer',...
      'Location','northEast');
% title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');ylabel('y-direction(km)');     % y축 이름

ti = 600;          % time interval 1hour = 3600sec
ub = u(1,:); vb = v(1,:);
stn = 2161; % neap 시기 유량 많을때 18일~20일
enn = 3168; 
sts = 3169; % spring 시기 유량 많을때 25일~27일
ens = 4176;
us_n = us(stn:enn); vs_n = vs(stn:enn);
um_n = um(stn:enn); vm_n = vm(stn:enn);
ub_n = ub(stn:enn); vb_n = vb(stn:enn);
us_s = us(sts:ens); vs_s = vs(sts:ens);
um_s = um(sts:ens); vm_s = vm(sts:ens);
ub_s = ub(sts:ens); vb_s = vb(sts:ens);

n = length(us_n);     % data의 size

% 시작점 위치 (0,0)
% 표층 data set 작성
nun_up = zeros(n,1); nvn_up =zeros(n,1);
for i=2:1:n
    nun_up(i) = nun_up(i-1) + us_n(i-1)/100000*ti;
    nvn_up(i) = nvn_up(i-1) + vs_n(i-1)/100000*ti;
end

% 중층 data set 작성
nun_mid = zeros(n,1); nvn_mid =zeros(n,1);
for i=2:1:n
    nun_mid(i) = nun_mid(i-1) + um_n(i-1)/100000*ti; 
    nvn_mid(i) = nvn_mid(i-1) + vm_n(i-1)/100000*ti;  
end
% 저층 data set 작성
nun_down = zeros(n,1); nvn_down =zeros(n,1);
for i=2:1:n
    nun_down(i) = nun_down(i-1) + ub_n(i-1)/100000*ti;
    nvn_down(i) = nvn_down(i-1) + vb_n(i-1)/100000*ti;
end
% 시작점 위치 (0,0)
% 표층 data set 작성
nus_up = zeros(n,1); nvs_up =zeros(n,1);
for i=2:1:n
    nus_up(i) = nus_up(i-1) + us_s(i-1)/100000*ti;
    nvs_up(i) = nvs_up(i-1) + vs_s(i-1)/100000*ti;
end

% 중층 data set 작성
nus_mid = zeros(n,1); nvs_mid =zeros(n,1);
for i=2:1:n
    nus_mid(i) = nus_mid(i-1) + um_s(i-1)/100000*ti; 
    nvs_mid(i) = nvs_mid(i-1) + vm_s(i-1)/100000*ti;  
end
% 저층 data set 작성
nus_down = zeros(n,1); nvs_down =zeros(n,1);
for i=2:1:n
    nus_down(i) = nus_down(i-1) + ub_s(i-1)/100000*ti;
    nvs_down(i) = nvs_down(i-1) + vb_s(i-1)/100000*ti;
end


subplot(143)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nun_up,nvn_up,'r');hold on;     % pgrogressive vector 도시
plot(nun_mid,nvn_mid,'g');hold on;
plot(nun_down,nvn_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
xlim([-100,100]);
ylim([-200,100]);
text(-70, 50,'Neap','fontsize',14);
axis equal
% legend('start','Upper layer','middle layer','bottom layer',...
%       'Location','northEast');
% title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');
% ylabel('y-direction(km)');     % y축 이름

subplot(144)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nus_up,nvs_up,'r');hold on;     % pgrogressive vector 도시
plot(nus_mid,nvs_mid,'g');hold on;
plot(nus_down,nvs_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
text(-70, 50,'Spring','fontsize',14);
xlim([-100,100]);
ylim([-200,100]);
axis equal
legend('start','Upper layer','middle layer','bottom layer',...
      'Location','northEast');
% title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');ylabel('y-direction(km)');     % y축 이름


%% progressive vector 표현 - 7일씩 neap/spring
% 유량이 적을때 각각 neap/spring
ti = 600;          % time interval 1hour = 3600sec
ub = u(1,:); vb = v(1,:);
stn = 2305; % neap 시기 유량 적을때 16일~18일
enn = 2592; 
sts = 3313; % spring 시기 유량 적을때
ens = 3600;
us_n = us(stn:enn); vs_n = vs(stn:enn);
um_n = um(stn:enn); vm_n = vm(stn:enn);
ub_n = ub(stn:enn); vb_n = vb(stn:enn);
us_s = us(sts:ens); vs_s = vs(sts:ens);
um_s = um(sts:ens); vm_s = vm(sts:ens);
ub_s = ub(sts:ens); vb_s = vb(sts:ens);
% 7일 간격 neap/spring
% stn = 145; % neap 시기 유량 적을때 7일간 2일~9일
% enn = 1152; 
% sts = 1153; % spring 시기 유량 적을때
% ens = 2160;
% us_n = us(stn:enn); vs_n = vs(stn:enn);
% um_n = um(stn:enn); vm_n = vm(stn:enn);
% ub_n = ub(stn:enn); vb_n = vb(stn:enn);
% us_s = us(sts:ens); vs_s = vs(sts:ens);
% um_s = um(sts:ens); vm_s = vm(sts:ens);
% ub_s = ub(sts:ens); vb_s = vb(sts:ens);
% n = length(us_n);     % data의 size

% 시작점 위치 (0,0)
% 표층 data set 작성
n = length(us_n);
nun_up = zeros(n,1); nvn_up =zeros(n,1);
for i=2:1:n
    nun_up(i) = nun_up(i-1) + us_n(i-1)/100000*ti;
    nvn_up(i) = nvn_up(i-1) + vs_n(i-1)/100000*ti;
end

% 중층 data set 작성
nun_mid = zeros(n,1); nvn_mid =zeros(n,1);
for i=2:1:n
    nun_mid(i) = nun_mid(i-1) + um_n(i-1)/100000*ti; 
    nvn_mid(i) = nvn_mid(i-1) + vm_n(i-1)/100000*ti;  
end
% 저층 data set 작성
nun_down = zeros(n,1); nvn_down =zeros(n,1);
for i=2:1:n
    nun_down(i) = nun_down(i-1) + ub_n(i-1)/100000*ti;
    nvn_down(i) = nvn_down(i-1) + vb_n(i-1)/100000*ti;
end
% 시작점 위치 (0,0)
% 표층 data set 작성
nus_up = zeros(n,1); nvs_up =zeros(n,1);
for i=2:1:n
    nus_up(i) = nus_up(i-1) + us_s(i-1)/100000*ti;
    nvs_up(i) = nvs_up(i-1) + vs_s(i-1)/100000*ti;
end

% 중층 data set 작성
nus_mid = zeros(n,1); nvs_mid =zeros(n,1);
for i=2:1:n
    nus_mid(i) = nus_mid(i-1) + um_s(i-1)/100000*ti; 
    nvs_mid(i) = nvs_mid(i-1) + vm_s(i-1)/100000*ti;  
end
% 저층 data set 작성
nus_down = zeros(n,1); nvs_down =zeros(n,1);
for i=2:1:n
    nus_down(i) = nus_down(i-1) + ub_s(i-1)/100000*ti;
    nvs_down(i) = nvs_down(i-1) + vb_s(i-1)/100000*ti;
end

figure('name','dry and wet')
set(gcf,'units', 'pixels', 'pos',[100 100 850 350])
subplot(141)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nun_up,nvn_up,'r');hold on;     % pgrogressive vector 도시
plot(nun_mid,nvn_mid,'g');hold on;
plot(nun_down,nvn_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
xlim([-100,100]);
ylim([-200,100]);
text(-70, 50,'Neap','fontsize',14);
axis equal
% legend('start','Upper layer','middle layer','bottom layer',...
%       'Location','northEast');
title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');ylabel('y-direction(km)');     % y축 이름

subplot(143)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nus_up,nvs_up,'r');hold on;     % pgrogressive vector 도시
plot(nus_mid,nvs_mid,'g');hold on;
plot(nus_down,nvs_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
text(-70, 50,'Spring','fontsize',14);
xlim([-100,100]);
ylim([-200,100]);
axis equal
legend('start','Upper layer','middle layer','bottom layer',...
      'Location','northEast');
% title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');ylabel('y-direction(km)');     % y축 이름

ti = 600;          % time interval 1hour = 3600sec
ub = u(1,:); vb = v(1,:);
stn = 2593; % neap 시기 유량 많을때 18일~20일
enn = 2880; 
sts = 3601; % spring 시기 유량 많을때 25일~27일
ens = 3888;
us_n = us(stn:enn); vs_n = vs(stn:enn);
um_n = um(stn:enn); vm_n = vm(stn:enn);
ub_n = ub(stn:enn); vb_n = vb(stn:enn);
us_s = us(sts:ens); vs_s = vs(sts:ens);
um_s = um(sts:ens); vm_s = vm(sts:ens);
ub_s = ub(sts:ens); vb_s = vb(sts:ens);

% 7일간격 neap(16일~23일) / spring(23일~30일)
% us_n = us(2161:3168); vs_n = vs(2161:3168);
% um_n = um(2161:3168); vm_n = vm(2161:3168);
% ub_n = ub(2161:3168); vb_n = vb(2161:3168);
% us_s = us(3888:4608); vs_s = vs(3888:4608);
% um_s = um(3888:4608); vm_s = vm(3888:4608);
% ub_s = ub(3888:4608); vb_s = vb(3888:4608);
% us_s = us(3169:4176); vs_s = vs(3169:4176);
% um_s = um(3169:4176); vm_s = vm(3169:4176);
% ub_s = ub(3169:4176); vb_s = vb(3169:4176);

n = length(us_n);     % data의 size

% 시작점 위치 (0,0)
% 표층 data set 작성
nun_up = zeros(n,1); nvn_up =zeros(n,1);
for i=2:1:n
    nun_up(i) = nun_up(i-1) + us_n(i-1)/100000*ti;
    nvn_up(i) = nvn_up(i-1) + vs_n(i-1)/100000*ti;
end

% 중층 data set 작성
nun_mid = zeros(n,1); nvn_mid =zeros(n,1);
for i=2:1:n
    nun_mid(i) = nun_mid(i-1) + um_n(i-1)/100000*ti; 
    nvn_mid(i) = nvn_mid(i-1) + vm_n(i-1)/100000*ti;  
end
% 저층 data set 작성
nun_down = zeros(n,1); nvn_down =zeros(n,1);
for i=2:1:n
    nun_down(i) = nun_down(i-1) + ub_n(i-1)/100000*ti;
    nvn_down(i) = nvn_down(i-1) + vb_n(i-1)/100000*ti;
end
% 시작점 위치 (0,0)
% 표층 data set 작성
nus_up = zeros(n,1); nvs_up =zeros(n,1);
for i=2:1:n
    nus_up(i) = nus_up(i-1) + us_s(i-1)/100000*ti;
    nvs_up(i) = nvs_up(i-1) + vs_s(i-1)/100000*ti;
end

% 중층 data set 작성
nus_mid = zeros(n,1); nvs_mid =zeros(n,1);
for i=2:1:n
    nus_mid(i) = nus_mid(i-1) + um_s(i-1)/100000*ti; 
    nvs_mid(i) = nvs_mid(i-1) + vm_s(i-1)/100000*ti;  
end
% 저층 data set 작성
nus_down = zeros(n,1); nvs_down =zeros(n,1);
for i=2:1:n
    nus_down(i) = nus_down(i-1) + ub_s(i-1)/100000*ti;
    nvs_down(i) = nvs_down(i-1) + vb_s(i-1)/100000*ti;
end


subplot(142)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nun_up,nvn_up,'r');hold on;     % pgrogressive vector 도시
plot(nun_mid,nvn_mid,'g');hold on;
plot(nun_down,nvn_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
xlim([-100,100]);
ylim([-200,100]);
text(-70, 50,'Neap','fontsize',14);
axis equal
% legend('start','Upper layer','middle layer','bottom layer',...
%       'Location','northEast');
% title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');
% ylabel('y-direction(km)');     % y축 이름

subplot(144)
% 시작점 기준위치 (0,0)
plot(0,0,'ks');hold on
plot(nus_up,nvs_up,'r');hold on;     % pgrogressive vector 도시
plot(nus_mid,nvs_mid,'g');hold on;
plot(nus_down,nvs_down,'b');hold on;
% 기준점에서 연결한 점선 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
text(-70, 50,'Spring','fontsize',14);
xlim([-100,100]);
ylim([-200,100]);
axis equal
legend('start','Upper layer','middle layer','bottom layer',...
      'Location','northEast');
% title('Progressive Vector Diagram at 3 Layer');          % 제목
% xlabel('x-direction(km)');ylabel('y-direction(km)');     % y축 이름




%% wave 의 종류 구분 standing/progressive
% x = [0:4745];
% depth = dd;
% figure()
% subplot(211)
% k = mu;
% [hAx,hLine1,hLine2]=plotyy(x,k/max(k), x,depth-mean(depth));
% % title('depth mean u velocity vs depth')
% ylabel(hAx(1),'u') % left y-axis
% ylabel(hAx(2),'depth') % right y-axis
% set(hAx(1),'xTick',[0:144:4746],'xlim',[0 4320],'xTickLabel',[0:1:33]);
% set(hAx(1),'ytick',[-3:3:3],'ylim',[-3 3]);
% set(hAx(2),'xTick',[0:144:4746],'xlim',[0 4320],'xTickLabel',[0:1:33]);
% set(hAx(2),'ytick',[-3:3:3],'ylim',[-3 3]);
% 
% subplot(212)
% s = mv;
% [hAx,hLine1,hLine2]=plotyy(x,s/max(s),x,depth-mean(depth));
% % title('depth mean v velocity vs depth')
% xlabel('Time (day)')
% ylabel(hAx(1), 'v') % left y-axis
% ylabel(hAx(2),'depth') % right y-axis
% set(hAx(1),'xTick',[0:144:4700],'xlim',[0 4320],'xTickLabel',[0:1:33]);
% set(hAx(1),'ytick',[-3:3:3],'ylim',[-3 3]);
% set(hAx(2),'xTick',[0:144:4700],'xlim',[0 4320],'xTickLabel',[0:1:33]);
% set(hAx(2),'ytick',[-3:3:3],'ylim',[-3 3]);
