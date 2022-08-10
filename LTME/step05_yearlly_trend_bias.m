
clc; clear all; close all

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
    %--------------------------------------------------

    figure('Position', [100, 100, 400, 400])
    subplot(211)
    hold on;
    plot(x1,yfit1,'k-.','linewidth',2);
    plot(x2,yfit2,'b-.','linewidth',2);
    plot(x1,y1,'k-','linewidth',2);
    plot(x2,y2,'b-','linewidth',2);
    
    legend('NFRDI','NOAA','location','north','Orientation','horizontal');
    legend('boxoff');
    text(8,-1.6,['Y_N_F = ' num2str(floor(P1(1)*1000)/1000) '*X + ' num2str(floor(P1(2)*1000)/1000) ]);
    text(8,-1,['Y_N_O = ' num2str(floor(P2(1)*1000)/1000) '*X + ' num2str(floor(P2(2)*1000)/1000) ]);
%     text(2,20,['Y_N2 = ' num2str(floor(P3(1)*10000)/10000) '*X + ' num2str(floor(P3(2)*1000)/1000) ]);
    title(num2str(fname(i)),'fontsize',14);
%     xlabel('Time (year)','fontsize',14);
    ylabel('Temp.','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 30]);
    set(gca,'ytick',[-2:2],'yticklabel',[-2:2],'ylim',[-2 2]);
    
    subplot(212)    
    er = (y2 - y1);  
    plot(x,er,'k','linewidth',2);
    line([1 30],[0 0]);
    title('Bias','fontsize',14);    
    xlabel('Time (year)','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 30]);
    ylabel('Temp.','fontsize',14);
    ylim([-2 2]);
    %--- save figure --------------------------------------
    outpath = '.\figure_yearly\';
    out=[outpath,num2str(fname(i)),'_bias'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=4;
    y_width=4;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
    %------------------------------------------------------
    
    close 
end



