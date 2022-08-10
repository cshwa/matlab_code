close all;clear all;clc


% 장기간의 yearly-mean SST trend 자료가 있는 경우 그림 그리기

%--- color map ------------------------------------------------------------
cmap = [255 85 0; 255 170 0;255 255 0;170 255 85;85 255 170; 0 255 255; ...
        0 170 255; 0 85 255;0 0 255];
cmap = flipud(cmap)/255;

%--- 사용할 자료 load -----------------------------------------------------
% load E:\11.사업\장기생태_3단계\draw_figure\data_from_yukyung\all_trend.txt
% data = all_trend(:,3);
%--- 아래 자료는 climate가 제외된 30년간의 daily mean 값이 있는 것이므로 trend
% 구하는 처리가 필요함
% load E:\11.사업\장기생태_3단계\draw_figure\data_from_yukyung\daily_trend_tongyeong.mat
% temp = daily_trend_tongyeong;
% y = nanmean(temp);
% x = 1:length(y);
% scatter(x,y);
% h = lsline;
% p1 = polyfit(get(h,'xdata'),get(h,'ydata'),1)

% reg_clim_no3
% reg_clim_temp
% reg_clim_salt

clearvars lon lat reg_clim_*
load('kodc_203_1997_2018_just_mean_yearly', 'coeff_*','lon_*','lat_*','reg_clim_*');
lon = lon_203; lat=lat_203; 
coef_s = coeff_salt; coef_t = coeff_temp; coef_n = coeff_no3; 
sp_no3_yr = reg_clim_no3; sp_temp_yr = reg_clim_temp; sp_salt_yr = reg_clim_salt;

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_* reg_clim_*
load('kodc_204_1997_2018_just_mean_yearly', 'coeff_*','lon_*','lat_*','reg_clim_*');
lon = cat(2,lon,lon_204); lat=cat(2,lat,lat_204);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3);
sp_no3_yr = cat(1,sp_no3_yr,reg_clim_no3); sp_temp_yr = cat(1,sp_temp_yr,reg_clim_temp); sp_salt_yr = cat(1,sp_salt_yr,reg_clim_salt);

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_* reg_clim_*
load('kodc_205_1997_2018_just_mean_yearly', 'coeff_*','lon_*','lat_*','reg_clim_*');
lon = cat(2,lon,lon_205); lat=cat(2,lat,lat_205);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3);
sp_no3_yr = cat(1,sp_no3_yr,reg_clim_no3); sp_temp_yr = cat(1,sp_temp_yr,reg_clim_temp); sp_salt_yr = cat(1,sp_salt_yr,reg_clim_salt);

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_* reg_clim_*
load('kodc_206_1997_2018_just_mean_yearly', 'coeff_*','lon_*','lat_*','reg_clim_*');
lon = cat(2,lon,lon_206); lat=cat(2,lat,lat_206);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3);
sp_no3_yr = cat(1,sp_no3_yr,reg_clim_no3); sp_temp_yr = cat(1,sp_temp_yr,reg_clim_temp); sp_salt_yr = cat(1,sp_salt_yr,reg_clim_salt);

clearvars coeff_salt coeff_temp coeff_no3 lon_* lat_* reg_clim_*
load('kodc_400_1997_2018_just_mean_yearly', 'coeff_*','lon_*','lat_*','reg_clim_*');
lon = cat(2,lon,lon_400); lat=cat(2,lat,lat_400);
coef_s = cat(1,coef_s,coeff_salt); coef_t =cat(1,coef_t, coeff_temp); coef_n = cat(1,coef_n,coeff_no3);
sp_no3_yr = cat(1,sp_no3_yr,reg_clim_no3); sp_temp_yr = cat(1,sp_temp_yr,reg_clim_temp); sp_salt_yr = cat(1,sp_salt_yr,reg_clim_salt);

return

%temp : - 0.04 ~ 0.04 
max(coef_t(:,1))
min(coef_t(:,1))

%salt : -0.03 ~ 0.03
max(coef_s(:,1))
min(coef_s(:,1))

%no3 : -1.6 ~ 0
max(coef_n(:,1))
min(coef_n(:,1))



%%temp
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
c_lim = 0.04;
clearvars data
data = coef_t(:,1);
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > 0.04);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','^oC/year')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SST trend (1997-2018)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_trend_temp_1997_2018'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------

%% salt
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
c_lim = 0.03;
clearvars data
data = coef_s(:,1);
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > 0.03);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','psu/year')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SSS trend (1997-2018)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_trend_salt_1997_2018'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------


%% no3
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
c_lim = 0;
clearvars data lon_plt
data = coef_n(:,1);
lon_plt = lon; lat_plt = lat;
lon_plt(data(:,1) == 0)=[];
lat_plt(data(:,1) == 0)=[];
data(data(:,1) == 0)=[];
for i = 1:length(cmap)-1
    c_int = [linspace(-1.6,0,9)];
    [a b] = find(data > 0);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','ug/L /year')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SSNO3 trend (1997-2018)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_trend_no3_1997_2018'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------


%% spatial_mean value & trends
nonan_case_n = find(isnan(sp_no3_yr(:,18))==0);
nonan_case_t = find(isnan(sp_temp_yr(:,18))==0);
nonan_case_s = find(isnan(sp_salt_yr(:,18))==0);
sp_mean_no3_yr = nanmean(sp_no3_yr(nonan_case_n,:),1);
sp_mean_temp_yr = nanmean(sp_temp_yr(nonan_case_t,:),1);
sp_mean_salt_yr = nanmean(sp_salt_yr(nonan_case_s,:),1);

% sp_mean_no3_yr
% sp_mean_temp_yr
% sp_mean_salt_yr

%% salt
clearvars yp_w_salt
color_pic = lines(size(sp_mean_salt_yr,1));
marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<',...
    'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<'};
xp = 1:40;
j=0
figure; hold on;
for i = 1:size(sp_mean_salt_yr,1)
  clearvars reg_data_salt xp_w_salt pf_w_salt
reg_data_salt = sp_mean_salt_yr(i,:);
% if isnan(reg_data_salt(1)) == 0 
    j = j+1;
    xp_w_salt = find(isnan(reg_data_salt)==0);
    pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
    yp_w_salt(i,:) = polyval(pf_w_salt,xp);
    scatter(1:40,sp_mean_salt_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:40, yp_w_salt(i,:),'color',color_pic(i,:));
    sp_coeff_salt(i,:) = pf_w_salt;
    sp_temp_case(j) = i;
end
hold on
% end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('정선관측-표층염분 연평균(공간평균)','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])


%% no3
clearvars yp_w_no3
j=0
figure; hold on;
for i = 1:size(sp_mean_no3_yr,1)
  clearvars reg_data_no3 xp_w_no3 pf_w_no3
reg_data_no3 = sp_mean_no3_yr(i,:);
% if isnan(reg_data_no3(1)) == 0 
    j = j+1;
    xp_w_no3 = find(isnan(reg_data_no3)==0);
    pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
    yp_w_no3(i,:) = polyval(pf_w_no3,xp);
    scatter(1:40,sp_mean_no3_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:40, yp_w_no3(i,:),'color',color_pic(i,:));
    sp_coeff_no3(i,:) = pf_w_no3;
    sp_temp_case(j) = i;
% end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('정선관측관측-표층질산염 연평균(공간평균)','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%temp
clearvars yp_w_temp
j=0
figure; hold on;
for i = 1:size(sp_mean_temp_yr,1)
  clearvars reg_data_temp xp_w_temp pf_w_temp
reg_data_temp = sp_mean_temp_yr(i,:);
% if isnan(reg_data_temp(1)) == 0 
    j = j+1;
    xp_w_temp = find(isnan(reg_data_temp)==0);
    pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
    yp_w_temp(i,:) = polyval(pf_w_temp,xp);
    scatter(1:40,sp_mean_temp_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:40, yp_w_temp(i,:),'color',color_pic(i,:));
    sp_coeff_temp(i,:) = pf_w_temp;
    sp_temp_case(j) = i;
% end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('정선관측관측-표층수온 연평균','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([13 20])



% absolute value increasing decreasing 

%%temp
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size

clearvars data
data = coef_t(:,1).*(length(1997:2018)-1); % rate to absolute value 
min(data)
max(data)
c_lim = 0.8;
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > c_lim);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','^oC')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SST 변화량 (1997-2018)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_abs_temp_1997_2018'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------

%% salt
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
clearvars data
data = coef_s(:,1).*(length(1997:2018)-1); % rate to absolute value ;
min(data)
max(data)
c_lim = 0.6;
for i = 1:length(cmap)-1    
    c_int = [-c_lim:(c_lim-(-c_lim))/8:c_lim];
    [a b] = find(data > c_lim);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon(a),lat(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','psu')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SSS 변화량 (1997-2018)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_abs_salt_1997_2018'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------


%% no3
%--- m_map 으로 coastal line load -----------------------------------------
figure1 = figure('Position', [100, 100, 500, 300]);
set(figure1,'PaperUnits','inches','PaperPosition',[0 0 1 1.5]*6);
m_proj('mercator','lon',[126.125 129.125],'lat',[33.125 35.125]);
m_gshhs_i('color','k');
m_gshhs_i('patch',[.9,.9,.9]);
m_grid('box','fancy','tickdir','in','linewidth',1);
hold on;

ms = 100; % marker size
clearvars data lon_plt
data = coef_n(:,1).*(length(1997:2018)-1); % rate to absolute value ;
lon_plt = lon; lat_plt = lat;
lon_plt(data(:,1) == 0)=[];
lat_plt(data(:,1) == 0)=[];
data(data(:,1) == 0)=[];
min(data)
max(data)
c_lim = 0;
for i = 1:length(cmap)-1
    c_int = [linspace(-35,0,9)];
    [a b] = find(data > c_lim);
%     m_plot(lon(a),lat(a),'MarkerEdgeColor',cmap(9,:), ...
%             'MarkerFaceColor',cmap(9,:),'MarkerSize',ms,'Marker','o');
    m_scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(9,:), ...
            'MarkerFaceColor',cmap(9,:),'Marker','o');
    [a b] = find(data > c_int(i) & data <= c_int(i+1));
%     m_plot(lon(a),lat(a),'Marker', 'o','MarkerEdgeColor',cmap(i,:), ...
%             'MarkerFaceColor',cmap(i,:),'MarkerSize',ms);
    m_scatter(lon_plt(a),lat_plt(a),ms,'MarkerEdgeColor',cmap(i,:), ...
            'MarkerFaceColor',cmap(i,:),'Marker', 'o');
hold on; 
end
cbh = colorbar;
set(cbh,'ytick',[0:1/8:1],'YTicklabel',c_int);
set(get(cbh,'Title'),'String','ug/L')
%--- 제목 수정하기 -------------------------------------
title('NFRDI SSNO3 변화량 (1997-2018)','fontsize',14);

%--- save figure --------------------------------------
    outpath = '.\figure\';
    %--- output 이름 수정하기 -------------------------
    out=['NFRDI_yearly_abs_no3_1997_2018'];
    set(gcf,'renderer','painter');
    set(gcf, 'PaperUnits', 'inches');
    x_width=5;
    y_width=3;
    set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
    saveas(gcf,out,'tif');
%------------------------------------------------------



