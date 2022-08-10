% 3������ ���� ���� ���Ӻ��͵��� �ۼ�
% input data: TC_data.dat
% m-file: progvec_exam2.m

clc;clear all;close all;
% Data �ҷ�����
[year,mm,dd,hr,min,u_up,v_up,u_mid,v_mid,...
   u_down,v_down]=textread('TC_data.dat','');
ti = 3600;          % time interval 1hour = 3600sec
n = length(year);     % data�� size

% ������ ��ġ (0,0)
% ǥ�� data set �ۼ�
nu_up = zeros(n,1); nv_up =zeros(n,1);
for i=2:1:n
    nu_up(i) = nu_up(i-1) + u_up(i-1)/100000*ti;
    nv_up(i) = nv_up(i-1) + v_up(i-1)/100000*ti;
end
% ���� data set �ۼ�
nu_mid = zeros(n,1); nv_mid =zeros(n,1);
for i=2:1:n
    nu_mid(i) = nu_mid(i-1) + u_mid(i-1)/100000*ti; 
    nv_mid(i) = nv_mid(i-1) + v_mid(i-1)/100000*ti;  
end
% ���� data set �ۼ�
nu_down = zeros(n,1); nv_down =zeros(n,1);
for i=2:1:n
    nu_down(i) = nu_down(i-1) + u_down(i-1)/100000*ti;
    nv_down(i) = nv_down(i-1) + v_down(i-1)/100000*ti;
end

% ������ ������ġ (0,0)
plot(0,0,'ks');hold on
plot(nu_up,nv_up,'r');hold on;     % pgrogressive vector ����
plot(nu_mid,nv_mid,'g');hold on;
plot(nu_down,nv_down,'b');hold on;

% ���������� ������ ���� 
plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-50 50],':k');
xlim([-400,50]);ylim([-450,50]);
axis equal
legend('start','Upper layer','middle layer','bottom layer',...
      'Location','SouthEast');
title('Progressive Vector Diagram at 3 Layer');          % ����
xlabel('x-direction(km)');ylabel('y-direction(km)');     % y�� �̸�