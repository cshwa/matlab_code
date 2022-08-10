%% progressive vector ǥ�� - ��ü�Ⱓ
% 3������ ���� ���� ���Ӻ��͵��� �ۼ�
% input data: TC_data.dat
% m-file: progvec_exam2.m

% clc;clear all;close all;
% Data �ҷ�����
% [year,mm,dd,hr,min,u_up,v_up,u_mid,v_mid,...
%    u_down,v_down]=textread('TC_data.dat','');
ti = 600;          % time interval 1hour = 3600sec
n = length(dd);     % data�� size
ub = u(1,:); vb = v(1,:);

stn = 1; % neap �ñ� ���� ������ 16��~18��
enn = en; 
% sts = 1153; % spring �ñ� ���� ������
% ens = 2160;

us_n = us(stn:enn); vs_n = vs(stn:enn);
um_n = um(stn:enn); vm_n = vm(stn:enn);
ub_n = ub(stn:enn); vb_n = vb(stn:enn);
% us_s = us(sts:ens); vs_s = vs(sts:ens);
% um_s = um(sts:ens); vm_s = vm(sts:ens);
% ub_s = ub(sts:ens); vb_s = vb(sts:ens);


% ������ ��ġ (0,0)
% ǥ�� data set �ۼ�
n = length(us_n);
nun_up = zeros(n,1); nvn_up =zeros(n,1);
for i=2:1:n
    nun_up(i) = nun_up(i-1) + us_n(i-1)/100000*ti;
    nvn_up(i) = nvn_up(i-1) + vs_n(i-1)/100000*ti;
end

% ���� data set �ۼ�
nun_mid = zeros(n,1); nvn_mid =zeros(n,1);
for i=2:1:n
    nun_mid(i) = nun_mid(i-1) + um_n(i-1)/100000*ti; 
    nvn_mid(i) = nvn_mid(i-1) + vm_n(i-1)/100000*ti;  
end
% ���� data set �ۼ�
nun_down = zeros(n,1); nvn_down =zeros(n,1);
for i=2:1:n
    nun_down(i) = nun_down(i-1) + ub_n(i-1)/100000*ti;
    nvn_down(i) = nvn_down(i-1) + vb_n(i-1)/100000*ti;
end

% ������ ��ġ (0,0)
% ǥ�� data set �ۼ�
% nus_up = zeros(n,1); nvs_up =zeros(n,1);
% for i=2:1:n
%     nus_up(i) = nus_up(i-1) + us_s(i-1)/100000*ti;
%     nvs_up(i) = nvs_up(i-1) + vs_s(i-1)/100000*ti;
% end
% 
% % ���� data set �ۼ�
% nus_mid = zeros(n,1); nvs_mid =zeros(n,1);
% for i=2:1:n
%     nus_mid(i) = nus_mid(i-1) + um_s(i-1)/100000*ti; 
%     nvs_mid(i) = nvs_mid(i-1) + vm_s(i-1)/100000*ti;  
% end
% % ���� data set �ۼ�
% nus_down = zeros(n,1); nvs_down =zeros(n,1);
% for i=2:1:n
%     nus_down(i) = nus_down(i-1) + ub_s(i-1)/100000*ti;
%     nvs_down(i) = nvs_down(i-1) + vb_s(i-1)/100000*ti;
% end

figure('name','weekly variation')
set(gcf,'units', 'pixels', 'pos',[100 100 350 350])
% subplot(141)
% ������ ������ġ (0,0)
plot(0,0,'ks');hold on
plot(nun_up,nvn_up,'r');hold on;     % pgrogressive vector ����
plot(nun_mid,nvn_mid,'g');hold on;
plot(nun_down,nvn_down,'b');hold on;
% ���������� ������ ���� 
plot([-250 250],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
plot([0 0],[-250 250],':k');
xlim([-150,150]);
ylim([-150,150]);
% text(-70, 50,'Neap','fontsize',14);
axis equal
legend('start','Upper layer','middle layer','bottom layer',...
      'Location','northEast');
title('Progressive Vector Diagram at 3 Layer');          % ����
% xlabel('x-direction(km)');ylabel('y-direction(km)');     % y�� �̸�

% subplot(142)
% % ������ ������ġ (0,0)
% plot(0,0,'ks');hold on
% plot(nus_up,nvs_up,'r');hold on;     % pgrogressive vector ����
% plot(nus_mid,nvs_mid,'g');hold on;
% plot(nus_down,nvs_down,'b');hold on;
% % ���������� ������ ���� 
% plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
% plot([0 0],[-50 50],':k');
% text(-70, 50,'Spring','fontsize',14);
% xlim([-100,100]);
% ylim([-200,100]);
% axis equal
% legend('start','Upper layer','middle layer','bottom layer',...
%       'Location','northEast');
% % title('Progressive Vector Diagram at 3 Layer');          % ����
% % xlabel('x-direction(km)');ylabel('y-direction(km)');     % y�� �̸�
% 
% ti = 600;          % time interval 1hour = 3600sec
% ub = u(1,:); vb = v(1,:);
% stn = 2161; % neap �ñ� ���� ������ 18��~20��
% enn = 3168; 
% sts = 3169; % spring �ñ� ���� ������ 25��~27��
% ens = 4176;
% us_n = us(stn:enn); vs_n = vs(stn:enn);
% um_n = um(stn:enn); vm_n = vm(stn:enn);
% ub_n = ub(stn:enn); vb_n = vb(stn:enn);
% us_s = us(sts:ens); vs_s = vs(sts:ens);
% um_s = um(sts:ens); vm_s = vm(sts:ens);
% ub_s = ub(sts:ens); vb_s = vb(sts:ens);
% 
% n = length(us_n);     % data�� size
% 
% % ������ ��ġ (0,0)
% % ǥ�� data set �ۼ�
% nun_up = zeros(n,1); nvn_up =zeros(n,1);
% for i=2:1:n
%     nun_up(i) = nun_up(i-1) + us_n(i-1)/100000*ti;
%     nvn_up(i) = nvn_up(i-1) + vs_n(i-1)/100000*ti;
% end
% 
% % ���� data set �ۼ�
% nun_mid = zeros(n,1); nvn_mid =zeros(n,1);
% for i=2:1:n
%     nun_mid(i) = nun_mid(i-1) + um_n(i-1)/100000*ti; 
%     nvn_mid(i) = nvn_mid(i-1) + vm_n(i-1)/100000*ti;  
% end
% % ���� data set �ۼ�
% nun_down = zeros(n,1); nvn_down =zeros(n,1);
% for i=2:1:n
%     nun_down(i) = nun_down(i-1) + ub_n(i-1)/100000*ti;
%     nvn_down(i) = nvn_down(i-1) + vb_n(i-1)/100000*ti;
% end
% % ������ ��ġ (0,0)
% % ǥ�� data set �ۼ�
% nus_up = zeros(n,1); nvs_up =zeros(n,1);
% for i=2:1:n
%     nus_up(i) = nus_up(i-1) + us_s(i-1)/100000*ti;
%     nvs_up(i) = nvs_up(i-1) + vs_s(i-1)/100000*ti;
% end
% 
% % ���� data set �ۼ�
% nus_mid = zeros(n,1); nvs_mid =zeros(n,1);
% for i=2:1:n
%     nus_mid(i) = nus_mid(i-1) + um_s(i-1)/100000*ti; 
%     nvs_mid(i) = nvs_mid(i-1) + vm_s(i-1)/100000*ti;  
% end
% % ���� data set �ۼ�
% nus_down = zeros(n,1); nvs_down =zeros(n,1);
% for i=2:1:n
%     nus_down(i) = nus_down(i-1) + ub_s(i-1)/100000*ti;
%     nvs_down(i) = nvs_down(i-1) + vb_s(i-1)/100000*ti;
% end
% 
% 
% subplot(143)
% % ������ ������ġ (0,0)
% plot(0,0,'ks');hold on
% plot(nun_up,nvn_up,'r');hold on;     % pgrogressive vector ����
% plot(nun_mid,nvn_mid,'g');hold on;
% plot(nun_down,nvn_down,'b');hold on;
% % ���������� ������ ���� 
% plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
% plot([0 0],[-50 50],':k');
% xlim([-100,100]);
% ylim([-200,100]);
% text(-70, 50,'Neap','fontsize',14);
% axis equal
% % legend('start','Upper layer','middle layer','bottom layer',...
% %       'Location','northEast');
% % title('Progressive Vector Diagram at 3 Layer');          % ����
% % xlabel('x-direction(km)');
% % ylabel('y-direction(km)');     % y�� �̸�
% 
% subplot(144)
% % ������ ������ġ (0,0)
% plot(0,0,'ks');hold on
% plot(nus_up,nvs_up,'r');hold on;     % pgrogressive vector ����
% plot(nus_mid,nvs_mid,'g');hold on;
% plot(nus_down,nvs_down,'b');hold on;
% % ���������� ������ ���� 
% plot([-50 50],[0 0],':k');   % plot([x1 x2],[y1 y2],'k') 
% plot([0 0],[-50 50],':k');
% text(-70, 50,'Spring','fontsize',14);
% xlim([-100,100]);
% ylim([-200,100]);
% axis equal
% legend('start','Upper layer','middle layer','bottom layer',...
%       'Location','northEast');
% % title('Progressive Vector Diagram at 3 Layer');          % ����
% % xlabel('x-direction(km)');ylabel('y-direction(km)');     % y�� �̸�


