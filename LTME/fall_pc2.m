% 2015 summer ����ķ ���� ���� �ڷ� ó�� 
% 2015.08.11.00:00 ~ 2015.08.27. 23:50

% clc; clear all; close all;
% ----------------------- Load xlsx data --------------------------------%
% directory = 'D:\mepl\data\01_Sumjin_yeosu_camp\2015_summer_ys\PC2_RCM9';
% filename = 'pc2.xlsx';
% filename = fullfile(directory,filename);
% data =  xlsread(filename);
% data = data(91:2555,:);
%-------------------------------------------------------------------------
% 
temp = readtable('fall_pc2.xlsx');
data = table2array(temp);
data = data(106:1866+106,:);
dir = data(:,2);
mag = data(:,1);


%%
% speed and dir to u and v
% mag = speed;
u = mag.*sind(dir);
v = mag.*cosd(dir);

%% depth �ڷ�� cutting
%{
for i = 1:length(depth)
    speed(i,depth(i)-1:end) = NaN;
    dir(i,depth(i)-1:end) = NaN;        
end

speed = flipud(rot90(speed(:,1:end)));
dir = flipud(rot90(dir(:,1:end)));

%}


%%
across = u;
up = v;
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 400, 300]);

subplot(211);
plot(across); 
% shading flat;
ylabel('Depth (m)','fontsize',14);
text(57,40, 'u (cm/s)');
ylim([-60 60]);
xlim([0 1866]);
title('Fall- PC2 Bottom current','fontsize',14);
set(gca, 'xtick',[1:144:13*144],'xticklabel',[22:30,1:31]);    
grid on;

subplot(212);
plot(up); 
% shading flat;
ylabel('Depth (m)','fontsize',14);
text(57,40, 'v (cm/s)');
% caxis([-80 80]);colorbar
set(gca, 'xtick',[1:144:13*144],'xticklabel',[22:30,1:31]);  
ylim([-60 60]);
xlim([0 1866]);
% set(gca,'yTick',[-1:2:15],'yTickLabel',[0:2:20],'ylim',[-1 14]); 
% hold on; plot(depth,'linewidth',2)
xlabel('Time (day)','fontsize',14);
grid on;



%% progressive vector
% Data �ҷ�����
u_vel = u;
v_vel = v;
ti = 600;          % time interval 1hour = 3600sec
n = length(u_vel);     % data�� size

% ������ ��ġ (0,0)
nu = zeros(n,1); nv =zeros(n,1);
for i=2:1:n
    nu(i) = nu(i-1) + u_vel(i-1)/100000*ti;   % nu�� u_vel�� ���� �� ����
    nv(i) = nv(i-1) + v_vel(i-1)/100000*ti;   % nv�� v_vel�� ���� �� ����
end
% ������ ������ġ (0,0)
  FigHandle = figure;
  set(FigHandle, 'Position', [200, 500,400, 400]);
plot(0,0,'ro');hold on
plot(nu,nv,'k');hold on;     % pgrogressive vector ����
plot([-50 50],[0 0],':b');   % plot([x1 x2],[y1 y2],'b') 
plot([0 0],[-50 50],':b');
xlim([-50,20]);
ylim([-20,50]);
axis equal
legend('start point','Bottom','Location','SouthEast');
title('Progressive Vector Diagram (Fall PC2)');   % ����
xlabel('x-direction(km)');ylabel('y-direction(km)');         % y�� �̸�

%%
% Scatterplot�� Compass grid�� ǥ���ϱ�
% m-file: compass_scatter.m

u_comp = u;
v_comp = v;
% ��� compas �׸���
figure()
polar(1,50);hold on  % compass�� ���ɿ��� ������ 100
% scatterplot ǥ��
plot(u_comp,v_comp,'b.','MarkerSize',5)   % marksize�� 5
set(gcf,'Color','w');       % ���� White


%%
% �ٴ��� ���� low pass filter
xb = u;
%--- low pass filter ��� �κ�---------------
del_t=1/6;
Nyq_fq=1/(2*del_t); 
cutoff_fq=1/25;  % 25: low pass filter ���� �ð�
wc=cutoff_fq/Nyq_fq;
order=8;
[b,a]=butter(order,wc,'low');
xb2=filtfilt(b,a,xb);
%--- low pass filter ��� �� ---------------------------------------------

FigHandle = figure;
set(FigHandle, 'Position', [100, 0, 500, 300]);
subplot(211)
plot(u);
hold on;
l = length(xb2);
plot(lpf_t:l-lpf_t,xb2(lpf_t:l-lpf_t),'k-','linewidth',1);
legend('u (cm/s)','Low pass filtered');
% plot(0:18*144,0,'b-');
set(gca,'xTick',[0:144:8928],'xTickLabel',[11:28,1:7],'xlim',[0 13*144]); 
set(gca,'ylim',[-60 60]);
grid on;
% xlabel('Time (day)','fontsize',14);
ylabel('cm/s','fontsize',14);
title('Fall PC2 Bottom current ','fontsize',20);

% �ٴ��� ���� low pass filter
xb =v;
%--- low pass filter ��� �κ�---------------
del_t=1/6;
Nyq_fq=1/(2*del_t); 
cutoff_fq=1/25;  % 25: low pass filter ���� �ð�
wc=cutoff_fq/Nyq_fq;
order=8;
[b,a]=butter(order,wc,'low');
xb2=filtfilt(b,a,xb);
%--- low pass filter ��� �� ---------------------------------------------

subplot(212)
plot(v);
l = length(xb2);
hold on;
plot(lpf_t:l-lpf_t,xb2(lpf_t:l-lpf_t),'k-','linewidth',1);
legend('u (cm/s)','Low pass filtered');
plot(xb2,'k:','linewidth',2);
% plot(0:18*144,0,'b-');
grid on;
set(gca,'xTick',[0:144:8928],'xTickLabel',[11:28,1:7],'xlim',[0 13*144]); 
set(gca,'ylim',[-60 60]);
xlabel('Time (day)','fontsize',14);
ylabel('cm/s','fontsize',14);
% title('Tidal residual current ','fontsize',20);

