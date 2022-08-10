clc;clear all;close all;

%--- color map ------------------------------------------------------------
% cmap = [0 0 255; 0 170 255; 255 255 0; 255 170 0; 255 85 0];
cmap = [122 16 228; 73 9 239; 10 92 253; 0 175 255; 139 249 255; 255 255 255; ...
        255 255 0; 255 192 0; 255 128 0; 255 65 0; 192 32 0; 128 0 0];
cmap = (cmap)/255;

%load file name
% name = dir('avhrr*');
ss = 1;
for mon = 2
mon=[num2char(mon,2)];
tn = 2014-1984+1;

    for yr = 1:tn
            filepath1 = 'E:\11.사업\장기생태_3단계\1차년도\Data\위성자료\monthly\';   %monthly folder
            name11 = ['avhrr_monthly' num2str(yr+1983) '_' mon '.nc'];
            file1 = strcat(filepath1, name11);
            nc = netcdf(file1); 
            temp1(yr,:,:) = nc{'temp'}(:);
            lon = nc{'long'}(:);
            lat = nc{'lat'}(:);
            
            temp2 = temp1(:,490:505,500:520);
            lon = lon(490:505,500:520);
            lat = lat(490:505,500:520);
    end
%     lim = [126 129.5 33 35.5];
%     [a b]=find(abs((lon(1,:)-lim(1)))==0.1250);
%     min_lon = min(b);    
%     [a b]=find(abs((lon(1,:)-lim(2)))==0.1250);
%     max_lon = max(b);    
%     [a b]=find(abs((lat(:,1)-lim(3)))<=0.1250);
%     min_lat = min(a);
%     [a b]=find(abs((lat(:,1)-lim(4)))<=0.1250);
%     max_lat = max(a);
    
%     temp2 = temp2(:,min_lat:max_lat,min_lon:max_lon);

    %trend data;
    x = [1:tn]';
    %--- edited by silver 원하는 영역만 계산하기 ---------
%     nu_lon = max_lon-min_lon+1;
%     nu_lat = max_lat-min_lat+1;
%     trend_data = zeros(nu_lat,nu_lon)+NaN;
    [a b]= size(lon);
    trend_data = zeros(a,b)+NaN;
    trend_idx = find(isnan(squeeze(temp2(1,:)))==0);
    %---------------------------------------------------
        for i=trend_idx
            data = temp2(:,i);
            [a,b] = polyfit(x,data,1);
            trend_data(i) = a(1);
        end
    [m n] = size(trend_data);
    T_trend_data(ss,1:m,1:n) = trend_data;

%%
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

m_pcolor(lon,lat,trend_data);
shading flat;
m_gshhs_i('patch',[.9,.9,.9]);
set(gca,'FontSize',14);
% tt = strcat('SST trend (1984-2014)');
% title(tt,'fontsize',20,'fontweight','bold');  
m_text(126.5,34.95,[num2str(mon) '월'],'color','k','FontSize',17,'fontweight','bold')
colormap(cmap);
caxis([-0.06,0.06]);
cbh = colorbar;
c_lim = 0.06;
c_min = -0.06;
c_int = [c_min:(c_lim-(c_min))/12:c_lim];
set(cbh,'ytick',c_int,'YTicklabel',c_int); 
ylabel(cbh,'C^{\circ}/ yr','fontsize',14);

ss = ss+1;
%--- save figure --------------------------------------
%     outpath = '.\figure\';
%     out=[outpath,'AVHRR_trend_',mon];
%     set(gcf,'renderer','painter');
%     set(gcf, 'PaperUnits', 'inches');
%     x_width=5;
%     y_width=3;
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
% ------------------------------------------------------

% clear temp2
% close all;
end
