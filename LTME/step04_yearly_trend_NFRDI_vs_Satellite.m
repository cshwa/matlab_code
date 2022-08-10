
clc; clear all; close all

%--- 정선별 정점 정보 가져오기 ---------------------------------------------
ST_info = load('position_silver.txt');
lat = ST_info(:,2);
lon = ST_info(:,3);
%--------------------------------------------------------------------------

%--- 정선별 정점 정보 가져오기 ---------------------------------------------
[ST_info] = textread('NSOposition_silver.txt','%s %*[^\n]','headerlines',1);
ST_info = str2mat(ST_info); 
ST_info = ST_info(:,[1:3,5:6]);
ST_info = str2num(ST_info);
fname = ST_info;
m = length(fname);

for i = 1:m
    fpath = ['./sorted_yearly/'];
    fname(i)
    sst_kodc = load([fpath, num2str(fname(i)),'_n_temp.dat']); % 수과원SST
    sst_sat = load([fpath, num2str(fname(i)),'_s_temp.dat']); % 위성SST
    
    %--- 자료들을 yearly mean 하기 ---------------------
    A = sst_kodc(:,2:end);
    [a b]= size(A);
    yearly_m_temp(i,1) = fname(i);
    if a < 30
        yearly_m_temp(i,2:a+1) = nanmean(A(1:a,:),2);
        sst_kodc = nanmean(A(1:a,:),2);
    else
        yearly_m_temp(i,2:31) = nanmean(A(1:30,:),2);
        sst_kodc = nanmean(A(1:30,:),2);
    end
      
    B = sst_sat(:,2:end);
    [a b]= size(B);
    yearly_m_temp(i,1) = fname(i);
    yearly_m_temp(i,2:31) = nanmean(B(1:30,:),2);
    sst_sat = nanmean(B(1:30,:),2);
    %---------------------------------------------------
    
    %--- linear fitting 구하기 -------------------------
    x = (1:30)';    
    y1 = sst_kodc;     x1 = x;    
    P1 = polyfit(x1,y1,1);
    yfit1 = P1(1)*x1+P1(2);
        
    y2 = sst_sat;     x2 = x;  
    P2 = polyfit(x2,y2,1);
    yfit2 = P2(1)*x2+P2(2);
    %--- trend 를 SST_trend 행렬로 정리하기 ------------
%     SST_trend(i,1) = fname(i);
%     SST_trend(i,2) = lat(i);
%     SST_trend(i,3) = lon(i);
    SST_t(i,1) = P1(1); % NFRDI trend
    SST_t(i,2) = P2(1); % NOAA/AVHRR trend
    %--------------------------------------------------

    figure('Position', [100, 100, 400, 200])
    hold on;
    plot(x1,yfit1,'k-','linewidth',2);
    plot(x2,yfit2,'b-','linewidth',2);
    plot(x1,y1,'k-.');
    plot(x2,y2,'b-.');
    
    legend('NFRDI','NOAA','location','north','Orientation','horizontal');
    legend('boxoff');
    text(8,-1.6,['Y_N_F = ' num2str(floor(P1(1)*1000)/1000) '*X + ' num2str(floor(P1(2)*1000)/1000) ]);
    text(8,-1,['Y_N_O = ' num2str(floor(P2(1)*1000)/1000) '*X + ' num2str(floor(P2(2)*1000)/1000) ]);
%     text(2,20,['Y_N2 = ' num2str(floor(P3(1)*10000)/10000) '*X + ' num2str(floor(P3(2)*1000)/1000) ]);
    title(num2str(fname(i)),'fontsize',14);
    xlabel('Time (year)','fontsize',14);
    ylabel('Temp.','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
    set(gca,'ytick',[-2:2],'yticklabel',[-2:2],'ylim',[-2 2]);
    
%     --- save figure --------------------------------------
%     outpath = '.\figure_yearly\';
%     out=[outpath,num2str(fname(i)),'yearly_trend'];
%     set(gcf,'renderer','painter');
%     set(gcf, 'PaperUnits', 'inches');
%     x_width=4;
%     y_width=2;
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
%     ------------------------------------------------------
    
    close 
end
NFRDI_t = SST_t(:,1);
AVHRR_t = SST_t(:,2);

save SST_trend.mat fname lon lat NFRDI_t AVHRR_t

