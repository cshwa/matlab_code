
% to_ye_temp 는 전체 영역에서 남해쪽만 뽑은 것임, 아래 정보 있음
% to_ye_temp(:,:,k) = data(488:502,501:520);

% kodc 정선과 일치하는 위성 자료
% 205-01,204-01,205-02 = 496,511
% 204-03,204-02 =496,510
% 204-05 = 495,509
% 204-04 = 496,509
% 204-06, 400-27 = 494,508
% 205-05 = 494,513
% 205-03 = 496,512
% 205-04 = 495,512
% 206-01 = 498,514
% 206-02,206-03 = 498,515
% 400-17 = 498,512
% 400-16,400-15 = 498,513
% 400-14 = 497,514
% 400-13 = 496,514
% 400-22 = 495,514
% 400-26 = 494,509
% 400-25 = 494,510
% 400-24 = 494,511

clc; clear all; 
% close all
noaa_r = [496 496 496 496 496 498 498 498 ...
            496 497 498 498 ...
            496 495 494 495 496 496 ...
            494 495 494 494 494];
noaa_c = [511 510 510 511 511 512 513 513 ...
            512 514 514 515 ...
            509 509 508 512 512 514 ...
            513 514 511 510 509];
kodc_n=[20401 20402 20403 20501 20502 40017 40016 40015 ...
        20503 40014 20601 20602 ...
        20404 20405 20406 20504 20603 40013 ...
        20505 40022 40024 40025 40026];
% 자료들을 연안과 외해로 구분    
noaa_r1 = noaa_r(1:12);
noaa_r2 = noaa_r(13:end);
noaa_c1 = noaa_c(1:12);
noaa_c2 = noaa_c(13:end);
kodc_n1 = kodc_n(1:12);
kodc_n2 = kodc_n(13:end);
%--- 위성자료 load --------------------------------------------------------
% 위성과  kodc 자료 84년부터 2015년까지 자료 사용하기
% 위성은 82년도 부터 있음 
load E:\11.사업\장기생태_3단계\Data\위성자료\from_YG\yearly\temp2;
sat_temp = temp2(1:31,:,:); % 84년부터 선택
load E:\11.사업\장기생태_3단계\Data\위성자료\from_YG\yearly\trend_data;
load E:\11.사업\장기생태_3단계\Data\위성자료\from_YG\yearly\trend_val;

%--- NFRDI 자료 load kodc(~1984) 자료 -------------------------------------
load E:\11.사업\장기생태_3단계\Data\국립수산과학원(KODC)_정선관측\yearly_m_temp
load E:\11.사업\장기생태_3단계\Data\국립수산과학원(KODC)_정선관측\yearly_temp

x = 1:31;
%--- kodc 자료에서 coastal/offshore -----------------------------
% inn = [1 2 3 7 8 9 12 13 16 17 18 19 ];
% ou = [4 5 6 10 11 14 15 20 21 22 23 24];
inn = [1 4 5 6 10 11 12 15 16 17 19 20 21 22];
ou = [2 3 7 8 9 13 14 18 23 24 25 26 27];
sst_kodc_inn = nanmean(yearly_m_temp(inn,2:end));
sst_kodc_ou = nanmean(yearly_m_temp(ou,2:end));
sst_sat_inn = nanmean(mean(sat_temp(:,noaa_r1,noaa_c1),2),3);
sst_sat_ou = nanmean(mean(sat_temp(:,noaa_r2,noaa_c2),2),3);
%----------------------------------------------------------------
%%
figure('Position', [100, 100, 900, 400])
%--- 연안 linear fitting -------
subplot(231)
    y1 = sst_kodc_inn;         
    P1 = polyfit(x,y1,1);
    yfit1 = P1(1)*x+P1(2);     
    
    y2 = sst_sat_inn;         
    P2 = polyfit(x,y2',1);
    yfit2 = P2(1)*x+P2(2);
    hold on;
    plot(x,y1,'k','linewidth',2);
    plot(x,y2,'b','linewidth',2);
    plot(x,yfit1,'k','linewidth',2);
    plot(x,yfit2,'b','linewidth',2);
    
    legend('NFRDI','NOAA','location','northwest','Orientation','horizontal');
    legend('boxoff');
    text(2,13.5,['Y_N_F = ' num2str(floor(P1(1)*1000)/1000) '*X + ' num2str(floor(P1(2)*100)/100) ]);
    text(2,15,['Y_N_O = ' num2str(floor(P2(1)*1000)/1000) '*X + ' num2str(floor(P2(2)*100)/100) ]);
    title('연안','fontsize',14);
%     xlabel('Time (year)','fontsize',14);
 ylabel('C^{\circ}','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
    set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
    

subplot(232)
    y3 = sst_kodc_ou;         
    P3 = polyfit(x,y3,1);
    yfit3 = P3(1)*x+P3(2);
       
    y4 = sst_sat_ou;         
    P4 = polyfit(x,y4',1);
    yfit4 = P4(1)*x+P4(2);
     hold on;
     
    plot(x,y3,'k','linewidth',2);
    plot(x,y4,'b','linewidth',2);
    plot(x,yfit3,'k','linewidth',2);
    plot(x,yfit4,'b','linewidth',2);
    
    text(2,13.5,['Y_N_F = ' num2str(floor(P3(1)*1000)/1000) '*X + ' num2str(floor(P3(2)*100)/100) ]);
    text(2,15,['Y_N_O = ' num2str(floor(P4(1)*1000)/1000) '*X + ' num2str(floor(P4(2)*100)/100) ]);
    title('외해','fontsize',14);
%     xlabel('Time (year)','fontsize',14);
%     ylabel('Temp.','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
    set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
    
subplot(233) % total 
    y5 = nanmean(yearly_m_temp(:,2:end));   
    P5 = polyfit(x,y5,1);
    yfit5 = P5(1)*x+P5(2);
       
    y6 = nanmean(mean(sat_temp(:,noaa_r,noaa_c),2),3);      
    P6 = polyfit(x,y6',1);
    yfit6 = P6(1)*x+P6(2);
     hold on;
     
    plot(x,y5,'k','linewidth',2);
    plot(x,y6,'b','linewidth',2);
    plot(x,yfit5,'k','linewidth',2);
    plot(x,yfit6,'b','linewidth',2);
    
    text(2,13.5,['Y_N_F = ' num2str(floor(P5(1)*1000)/1000) '*X + ' num2str(floor(P5(2)*100)/100) ]);
    text(2,15,['Y_N_O = ' num2str(floor(P6(1)*1000)/1000) '*X + ' num2str(floor(P6(2)*100)/100) ]);
    title('남해 전체','fontsize',14);
%     xlabel('Time (year)','fontsize',14);
%     ylabel('Temp.','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
    set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
    
    %--- bias 계산 ---------------
    subplot(234)  
    er = y2'-y1;
    P = polyfit(x,er,1);
    yfit = P(1)*x+P(2);
    plot(x,er,'k','linewidth',2);    
    xlabel('Time (year)','fontsize',14);
     ylabel('C^{\circ}','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
    set(gca,'ytick',[-2:2:6],'yticklabel',[-2:2:6],'ylim',[-2 3]);
%     text(2, -1, ['Mean: ' num2str(floor(mean(y2'-y1)*100)/100)],'fontsize',14);
    grid on;
    
    subplot(235) 
    er = y4'-y3;
    P = polyfit(x,er,1);
    yfit = P(1)*x+P(2);
    plot(x,er,'k','linewidth',2);
    xlabel('Time (year)','fontsize',14);
%     ylabel('Temp.','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
    set(gca,'ytick',[-2:2:6],'yticklabel',[-2:2:6],'ylim',[-2 3]);
%     text(2, -1, ['Mean: ' num2str(floor(mean(y4'-y3)*100)/100)],'fontsize',14);
    grid on; 
    
    subplot(236)
    er = y6'-y5;
    P = polyfit(x,er,1);
    yfit = P(1)*x+P(2);
%     plot(x,er,'k',x,yfit,'k','linewidth',2);
     plot(x,er,'k','linewidth',2);
    xlabel('Time (year)','fontsize',14);
%     ylabel('Temp.','fontsize',14);
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
    set(gca,'ytick',[-2:2:6],'yticklabel',[-2:2:6],'ylim',[-2 3]);
%     text(2, -1, ['Mean: ' num2str(floor(mean(y6'-y5)*100)/100)],'fontsize',14);
    grid on; 

    %--- save figure --------------------------------------
    outpath = '.\figure_yearly\';
    out=[outpath,'SST_inn_ou'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=9;
    y_width=4;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
    %------------------------------------------------------



