%===================================================================
% ����������п�(KODC) �������� �ڷ� - ���縸 �߽� �ڷ� ����ȭ (����)
% 
% 
% �� �ܰ���� ó���� ���� �ڷ�� ���ϱ� ���� �����ڷ� ������ �Ѿ �۾���
% copywrite EunByeol Cho 2015.05.24
%===================================================================

clc; clear all; close all

load yearly_m_temp_o_s;
data = yearly_m_temp_o_s;
[a b]= size(data);
for i = 1:a
    y = data(i,2:32); %1984-2014����� yearly-mean SST
    x = 1:31;
    p(i,:) = polyfit(x,y,1);
    NFRDI_t(i,1) = p(1);
    scatter(x,y);
    h = lsline;
    p1(i,:) = polyfit(get(h,'xdata'),get(h,'ydata'),1)
end


%--- ������ ���� ���� �������� ---------------------------------------------
[ST_info] = textread('NSOposition_silver.txt','%s %*[^\n]','headerlines',1);
ST_info = str2mat(ST_info); 
ST_info = ST_info(:,[1:3,5:6]);
ST_info = str2num(ST_info);
fname = ST_info;

m = length(fname);
yearly_temp = zeros(m,6*31+1);
for ss = 1:length(ST_info)
% for ss = 1:1
    fpath = ['.\sorted_yearly\'];
    data = load([fpath, num2str(fname(ss)),'_1984_o_s.dat']);
    A = data(:,2:end);
    [a b]= size(A);
    
    % ¦���� ������ 1984���� 2014���� ǥ������
    B = reshape(A',[a*b,1]);
    yearly_temp_o_s(ss,1) = fname(ss);
    yearly_temp_o_s(ss,2:(6*31+1)) = B(1:6*31);
    % �� ��� ������ trend��
    yearly_m_temp_o_s(ss,1) = fname(ss);
    yearly_m_temp_o_s(ss,2:32) = nanmean(A(1:31,:),2);
end



%% �� ��� ���� trend
figure('Position', [100, 100, 700, 300])
plot(1:31,yearly_m_temp(:,2:end)');
hold on;
scatter(1:31,nanmean(yearly_m_temp(:,2:end)),'linewidth',2);
h=lsline;
set(h,'LineWidth',2,'color','k')
xlabel('Time (year)','fontsize',14);
ylabel('Temp.','fontsize',14);
set(gca,'xtick',[1:5:30],'xticklabel',[1984:5:2015],'xlim',[1 31]);
% set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
% legend('20401','20402','20403','20404','location','northeast');
grid on;
title('South sea - coastal SST','fontsize',14);
%--- save figure --------------------------------------
% outpath = '.\figure\';
% out=[outpath,'yearly_trend_from1984'];
% set(gcf,'renderer','painter');
% set(gcf, 'PaperUnits', 'inches');
% x_width=7;
% y_width=3;
% set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
% saveas(gcf,out,'tif');
%------------------------------------------------------

%% �� ��� ���� trend
figure('Position', [100, 100, 700, 300])
for i = 1:m
scatter(1:31,yearly_m_temp(i,2:end)');
hold on;
lsline
end
xlabel('Time (year)','fontsize',14);
ylabel('Temp.','fontsize',14);
set(gca,'xtick',[1:5:30],'xticklabel',[1985:5:2015],'xlim',[1 33]);
% set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
% legend('20401','20402','20403','20404','location','northeast');
grid on;