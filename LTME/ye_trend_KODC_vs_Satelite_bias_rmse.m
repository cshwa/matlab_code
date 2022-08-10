
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
noaa_r = [496 496 496 496 496 498 498 498 ...
            496 497 498 498 ...
            496 495 494 495 496 496 ...
            494 495 494 494 494 494];
noaa_c = [511 510 510 511 511 512 513 513 ...
            512 514 514 515 ...
            509 509 508 512 512 514 ...
            513 514 511 510 509 508];
kodc_n=[20401 20402 20403 20501 20502 40017 40016 40015 ...
        20503 40014 20601 20602 ...
        20404 20405 20406 20504 20603 40013 ...
        20505 40022 40024 40025 40026 40027];

% 위성과  kodc 자료 84년부터 2015년까지 자료 사용하기
% 위성은 82년도 부터 있음 
load to_ye_temp
sat_temp = to_ye_temp(:,:,3:end-2); % 84년부터 선택

% load kodc(~1984) 자료 
load E:\11.사업\장기생태_3단계\Data\국립수산과학원(KODC)_정선관측\yearly_m_temp
% load E:\11.사업\장기생태_3단계\Data\국립수산과학원(KODC)_정선관측\yearly_temp

%--- 한 그림에 여러개의 정점 ------------------
figure('Position', [100, 100, 1500, 400])
%--------------------------------------------
x = 1:31;
k = 1;
l =length(noaa_r);
for i = 1:l
    [a b]= find(yearly_m_temp(:,1) == kodc_n(i));
    sst_kodc = yearly_m_temp(a,2:end);
    sst_sat = squeeze(sat_temp(noaa_r(i)-488,noaa_c(i)-501,:));
        
    y1 = sst_kodc;   x1 = x;      
    y2 = sst_sat';   x2 = x;    
   
    er = (y2 - y1);    % Errors
    er_all(i,:) = er;
    ser = (y2 - y1).^2;   % Squared Error
    mser = mean((y2 - y1).^2);   % Mean Squared Error
    RMSE = sqrt(mean((y2 - y1).^2));  % Root Mean Squared Error
    
    h=subplot(3,8, k);
    p = get(h, 'pos');
    p(3) = p(3) + 0.01; %[left, bottom, width, height]
    p(4) = p(4) + 0.015;
    set(h, 'pos', p);
%     figure('Position', [100, 100, 300, 200])
    hold on;

%     scatter(x1,er,'k.');
    %--- linear fitting --------------
    y1 = er;            
    P = polyfit(x1,y1,1);
    yfit = P(1)*x1+P(2);
    plot(x1,er,'b',x1,yfit,'b','linewidth',2);
%     h=lsline;
%     set(h,'linewidth',2,'color','k');
    line([1 31],[0 0]);
%     grid on;

    legend('boxoff');
    text(3, 4, [num2str(kodc_n(i))],'fontsize',12);
%     text(3, 0.8, [num2str((round(P(1)*10000))/10000)],'fontsize',12);
    text(3, -1.2, [num2str((round(P(1)*1000))/1000) '/' num2str((round(P(2)*100))/100)],'fontsize',12);
    if k > 16
        xlabel('Time (year)','fontsize',14);
    end
    if mod(k,8)==1
        ylabel('Temp.','fontsize',14);
    end
    set(gca,'xtick',[1:8:30],'xticklabel',[1984:8:2015],'xlim',[1 29]);
    set(gca,'ytick',[-2:2:6],'yticklabel',[-2:2:6],'ylim',[-2 5]);
    
    %--- 한 정점씩 save figure -----------------------------
%     outpath = '.\figure_yearly\';
%     out=[outpath,num2str(kodc_n(i)),'_bias'];
%     set(gcf,'renderer','painter');
%     set(gcf, 'PaperUnits', 'inches');
%     x_width=3;
%     y_width=2;
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
    %------------------------------------------------------
%     close 
    
%     close 
    k = k+1;
end

%--- save figure --------------------------------------
    outpath = '.\figure_yearly\';
    out=[outpath,'all_bias'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=15;
    y_width=4;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------

