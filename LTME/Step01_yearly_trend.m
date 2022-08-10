%===================================================================
% 국립수산과학원(KODC) 정선관측 자료 - 광양만 중심 자료 가시화 (남해)
% 
% 
% 
% copywrite EunByeol Cho 2015.05.24
%===================================================================

clc; clear all; close all

%--- 정선별 정점 정보 가져오기 ---------------------------------------------
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
    data = load([fpath, num2str(fname(ss)),'_1984.dat']);
    A = data(:,2:end);
    [a b]= size(A);
    
    % 짝수달 순서로 1984부터 2014까지 표층수온
    B = reshape(A',[a*b,1]);
    yearly_temp(ss,1) = fname(ss);
    yearly_temp(ss,2:(6*31+1)) = B(1:6*31);
    % 년 평균 수온을 trend로
    yearly_m_temp(ss,1) = fname(ss);
    yearly_m_temp(ss,2:32) = nanmean(A(1:31,:),2);
end

save yearly_m_temp yearly_m_temp
save yearly_temp yearly_temp
%% 년 평균 수온 trend
figure('Position', [100, 100, 700, 300])
plot(1:31,yearly_m_temp(:,2:end)');
hold on;
scatter(1:31,nanmean(yearly_m_temp(:,2:end)),'linewidth',2);
h=lsline;
set(h,'LineWidth',2,'color','k')
xlabel('Time (year)','fontsize',14);
ylabel('Temp.','fontsize',14);
set(gca,'xtick',[1:5:30],'xticklabel',[1984:5:2015],'xlim',[1 31]);
set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
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

%% 년 평균 수온 trend
%{
figure('Position', [100, 100, 700, 300])
for i = 1:22
scatter(1:31,yearly_m_temp(i,2:end)');
hold on;
lsline
end
xlabel('Time (year)','fontsize',14);
ylabel('Temp.','fontsize',14);
set(gca,'xtick',[1:5:30],'xticklabel',[1985:5:2015],'xlim',[1 33]);
set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
legend('20401','20402','20403','20404','location','northeast');
grid on;
%}