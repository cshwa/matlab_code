
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

clc; clear all; close all
% noaa_r=[496 496 496 496 496 495 496 494 494 494 496 495 498 498 498 498 498 498 497 496 495 494 494 494];
% noaa_c=[511 511 510 511 510 509 509 508 508 513 512 512 514 515 515 512 513 513 514 514 514 509 510 511];
% kodc_n=[20501 20401 20402 20502 20403 20405 20404 20406 40027 20505 20503 20504 20601 ...
%     20602 20603 40017 40016 40015 40014 40013 40022 40026 40025 40024];
noaa_r = [496 495 495 494 496 496 496 496 496 498 498 498 ...
            496 497 498 498 ...
            496 495 494 495 496 496 ...
            494 495 494 494 494];
noaa_c = [506 506 506 509 511 510 510 511 511 512 513 513 ...
            512 514 514 515 ...
            509 509 508 512 512 514 ...
            513 514 511 510 509];
kodc_n=[20301 20302 20303 40027 20401 20402 20403 20501 20502 40017 40016 40015 ...
        20503 40014 20601 20602 ...
        20404 20405 20406 20504 20603 40013 ...
        20505 40022 40024 40025 40026];

% 위성과  kodc 자료 84년부터 2015년까지 자료 사용하기
% 위성은 82년도 부터 있음 
load E:\11.사업\장기생태_3단계\Data\위성자료\from_YG\yearly\temp2;
sat_temp = temp2(1:31,:,:); % 84년부터 선택
load E:\11.사업\장기생태_3단계\Data\위성자료\from_YG\yearly\trend_data;
load E:\11.사업\장기생태_3단계\Data\위성자료\from_YG\yearly\trend_val;

% load kodc(~1984) 자료 
load E:\11.사업\장기생태_3단계\Data\국립수산과학원(KODC)_정선관측\yearly_m_temp
load E:\11.사업\장기생태_3단계\Data\국립수산과학원(KODC)_정선관측\yearly_temp

x = 1:31;
l =length(noaa_r);
for i = 1:4
    [a b]= find(yearly_m_temp(:,1) == kodc_n(i));
    sst_kodc = yearly_m_temp(a,2:end);
    sst_sat = squeeze(sat_temp(:,noaa_r(i),noaa_c(i)));
        
%     y1 = sst_kodc;     x1 = x;    
%     P1 = polyfit(x1,y1,1);
%     yfit1 = P1(1)*x1+P1(2);
    y2 = sst_sat';     x2 = x;  
    P2(1) = trend_data(noaa_r(i),noaa_c(i));
    P2(2) = trend_val(noaa_r(i),noaa_c(i));
    yfit2 = P2(1)*x2+P2(2);
    noaa(i,1) = P2(1);
    %--- linear fitting --------------
%     P1 = polyfit(x1,y1,1);
%     yfit1 = P1(1)*x1+P1(2);
%     P3 = polyfit(x2,y2,1);
%     yfit3 = P3(1)*x2+P3(2);

%     figure('Position', [100, 100, 300, 200])
%     hold on;
%     plot(x1,yfit1,'k-','linewidth',2);
%     plot(x2,yfit2,'b-','linewidth',2);
% %     plot(x2,yfit3,'g-','linewidth',2);
%     plot(x1,y1,'k-.');
%     plot(x2,y2,'b-.');
%     legend('KODC','NOAA','location','north','Orientation','horizontal');
%     legend('boxoff');
%     text(2,13.5,['Y_K = ' num2str(floor(P1(1)*1000)/1000) '*X + ' num2str(floor(P1(2)*1000)/1000) ]);
%     text(2,15,['Y_N = ' num2str(floor(P2(1)*1000)/1000) '*X + ' num2str(floor(P2(2)*1000)/1000) ]);
% %     text(2,20,['Y_N2 = ' num2str(floor(P3(1)*10000)/10000) '*X + ' num2str(floor(P3(2)*1000)/1000) ]);
%     title(num2str(kodc_n(i)),'fontsize',14);
%     xlabel('Time (year)','fontsize',14);
%     ylabel('Temp.','fontsize',14);
%     set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 31]);
%     set(gca,'ytick',[13:2:23],'yticklabel',[13:2:23],'ylim',[13 23]);
    %--- save figure --------------------------------------
%     outpath = '.\figure_yearly\';
%     out=[outpath,num2str(kodc_n(i)),'_o_s'];
%     set(gcf,'renderer','painter');
%     set(gcf, 'PaperUnits', 'inches');
%     x_width=3;
%     y_width=2;
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
    %------------------------------------------------------
    close 
end



