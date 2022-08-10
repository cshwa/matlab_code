%% Satillite chl
close all; clear; clc;

% GMIS_S_CHLA_01_1998.nc
var_name = 'chl';
time_dyear = 1997:2005;

% f_name = 'hawaii_soest_2002-2005_v5.1.nc';  % satellite product SEAWIF seawif monthly v5.1
f_name = 'hawaii_soest_1997-2001.nc'; % satellite product SEAWIF seawif monthly v5.2
f_name2 = 'hawaii_soest_2002-2005.nc'; % satellite product SEAWIF seawif monthly v5.2
chl_1=ncread(f_name,'ch'); % 1997-09-16T00:00:00Z ~ 2002-01-16T00:00:00Z
chl_2=ncread(f_name2,'ch'); % 1997-09-16T00:00:00Z ~ 2002-01-16T00:00:00Z
lon=ncread(f_name,'longitude');
lat=ncread(f_name,'latitude');
t1=ncread(f_name,'time');
t2=ncread(f_name2,'time');

[y x]=meshgrid(lat,lon);

area_mask=NaN(size(x,1),size(x,2));
area_mask(find(x >= 126.5 & x <= 130 & 32.5<= y & 35.4 >= y))=1;
% area_mask(index_nan) = NaN;

% chl_1(:,:,end)=[];
chl_2(:,:,end)=[];
chl_d = [chl_1];
chl_d(:,:,end+1:end+size(chl_2,3)) = chl_2;

t_d = [t1];
t_d(end+1:end+length(t2)-1) = t2(1:end-1);

% pcolor(x,y,squeeze(nanmean(chl_d,3))); shading flat;
figure;
pcolor(x,y,squeeze(chl_d(:,:,1))); shading flat;

t_d_s=utime2date(t_d);
datestr(t_d_s)

d_mean_pre = chl_d(:,:,5:end);

for i = 1:12
    d_mean_pre2(:,:,i)=nanmean(d_mean_pre(:,:,i:12:end),3);
end

d_mean(:,:,1:8) = d_mean_pre2(:,:,1:8);
k = 0;
for i = 9:12
    k=k+1;
d_mean(:,:,i)=nanmean(cat(3,chl_d(:,:,k),d_mean_pre2(:,:,i)),3);
end

for i= 1:12
spatial_mean(i,1) = nanmean(nanmean(squeeze(d_mean(:,:,i)) .* (area_mask./area_mask),1),2);
end

for i= 1:4
spatial_season(i,1) = nanmean(nanmean(squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) .* (area_mask./area_mask),1),2);
end

save('seaWIFS_spmean_monthly.mat','spatial_mean','spatial_season');

for i = 1:12
    fig=figure;   
%         pcolor(x,y,squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) ); shading flat; 
        pcolor(x,y,squeeze(d_mean(:,:,i))); shading flat; 
        h=colorbar; caxis([0 10]);
    title([var_name ,' 0.25  seaWIFS v5.2 monthly ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' num2str(i), '-mon']);
    ylim([15 52.0063])
    xlim([115 161.90])
    plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
    ylim([15 52.0063])
    xlim([115 161.90])
    grid on;
    colormap('jet'); set(get(h,'title'),'string','Chl.a (ug/L)');
    saveas(fig,[var_name ,'_0.25_seaWIFS_v5.2_monthly_mean_', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' num2str(i), '-mon.png']);
    close;
end

seasonal = {'winter','spring','summer','fall'};
for i = 1:4
    fig=figure;   
        pcolor(x,y,squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) ); shading flat; 
%         pcolor(x,y,squeeze(d_mean(:,:,i))); shading flat; 
        h=colorbar; caxis([0 10]);
    title([var_name ,' 0.25  seaWIFS v5.2 monthly ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' seasonal{i}]);
    ylim([15 52.0063])
    xlim([115 161.90])
    plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
    ylim([15 52.0063])
    xlim([115 161.90])
    grid on;
    colormap('jet'); set(get(h,'title'),'string','Chl.a (ug/L)');
    saveas(fig,[var_name ,'_0.25_seaWIFS_v5.2_seasonal_', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' seasonal{i}, '.png']);
    close;
end


nan_tgt=squeeze(nanmean(d_mean,3));
save('seaWIFS_v5.2_nantgt.mat','nan_tgt');

% close all; clear; clc;
% 
% var_name = 'chl';
% chl_1=ncread('erdSW2018chlamday_1997-2001.nc','chlorophyll'); % 1997-09-16T00:00:00Z ~ 2002-01-16T00:00:00Z
% chl_2=ncread('erdSW2018chlamday_2002-2005.nc','chlorophyll'); % 2002-01-16T00:00:00Z ~ 2006-01-16T00:00:00Z
% lon=ncread('erdSW2018chlamday_1997-2001.nc','longitude');
% lat=ncread('erdSW2018chlamday_1997-2001.nc','latitude');
% t_d_1=ncread('erdSW2018chlamday_1997-2001.nc','time');
% t_d_2=ncread('erdSW2018chlamday_2002-2005.nc','time');
% 
% [y x]=meshgrid(lat,lon);
% 
% t_d_1(end)=[]; t_d_2(end)=[]; 
% chl_1(:,:,end)=[]; chl_2(:,:,end)=[];
% chl_d = [chl_1];
% chl_d(:,:,end+1:end+size(chl_2,3)) = chl_2;
% 
% t_d = [t_d_1];
% t_d(end+1:end+length(t_d_2)) = t_d_2;
% 
% % pcolor(x,y,squeeze(nanmean(chl_d,3))); shading flat;
% figure;
% pcolor(x,y,squeeze(chl_d(:,:,1))); shading flat;
% 
% t_d_s=utime2date(t_d);
% datestr(t_d_s)
% 
% 
% d_mean_pre = chl_d(:,:,5:end);
% 
% 
% 
% if plt_on == 1
%     fig=figure;   
%         pcolor(lon,lat,squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) ); shading flat; 
%         colorbar; caxis([0 10]);
%     title([var_name ,' interp ' mod_name ' historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' num2str(i), '-mon']);
%     ylim([15 52.0063])
%     xlim([115 161.90])
%     plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
%     ylim([15 52.0063])
%     xlim([115 161.90])
%     colormap('jet')
%     saveas(fig,[var_name ,' interp ' mod_name ' historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' num2str(i), '-mon.png']);
%     close;
% end



close all; clear; clc;
% monthly_mean(extract_other_seas) = NaN;

% time_dyear = 1965:2005; 
time_dyear = 1997:2005; 
var_name = 'chl'
% var_name = 'no3'
% var_name = 'o2'
mod_name_pre = {'IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR','GFDL-ESM2M', ...
    'GFDL-ESM2G','CNRM-CM5','CESM1-BGC'}; 
% mod_name_pre = {'CNRM-CM5'}; 
plt_on = 0; %switch to plot smth.

% standard_name = 'mole_concentration_of_nitrate_in_sea_water'
% long_name     = 'Dissolved Nitrate Concentration at Surface'
% units         = 'mol m-3'
% 1 mol/m3 = 1000 mmol/m3

x=ncread([var_name '_merged_CESM1_BGC_his_1965-2005.nc'],'lon');
y=ncread([var_name '_merged_CESM1_BGC_his_1965-2005.nc'],'lat');

[lat lon]=meshgrid(y,x);

area_mask=NaN(size(lon,1),size(lon,2));
area_mask(find(lon >= 126.5 & lon <= 130 & 32.5<= lat & 35.4 >= lat))=1;
% area_mask(index_nan) = NaN;
extract_other_seas=find(isnan(area_mask)==1);
seasonal = {'winter','spring','summer','fall'};

for l = 1:length(mod_name_pre)
    clearvars -except time_dyear var_name mod_name_pre l lon lat spatial_mean extract_other_seas plt_on area_mask seasonal spatial_season

mod_name = mod_name_pre{l};
disp(mod_name)
%     D:\장기생태\Dynamic\CMIP5_bio\no3
cd(['D:\장기생태\Dynamic\CMIP5_bio\' var_name])

load([var_name ,'_interp_' mod_name '_historical_clim_', num2str(time_dyear(1)), '-', num2str(time_dyear(end)), '.mat'],'d_mean');

for i = 1:4
    
%     d_mean(d_mean == 0) = NaN;
if plt_on == 1
    fig=figure;    
    switch(var_name)
        case 'chl'
        pcolor(lon,lat,squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) .*10^6); %kg/m^3 to ug/L 
        shading flat; 
        h=colorbar; caxis([0 10]);
        colorbarlabel='Chl (ug/L)';
        case 'no3'
        pcolor(lon,lat,squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) .* 1000); shading flat; 
        h=colorbar; caxis([0 10]);
        colorbarlabel='no_3 (mmol/m^3)';
        case 'o2'
        pcolor(lon,lat,squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) .* 1000); shading flat;     
        h=colorbar; caxis([0 380]);
        colorbarlabel='O_2 (mmol/m^3)';
    end
    title([var_name ,' interp ' mod_name ' historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' seasonal{i}]);
    ylim([15 52.0063])
    xlim([115 161.90])
    plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
    ylim([15 52.0063])
    xlim([115 161.90])
    colormap('jet'); set(get(h,'title'),'string',colorbarlabel);
    saveas(fig,[var_name ,' interp ' mod_name ' historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' seasonal{i} '.png']);
    close;
end   
%%  spatial mean
%     fig=figure;    
    switch(var_name)
        case 'chl'
            clearvars temp_d
        temp_d = (squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) .* 10^6) ; shading flat; 
        spatial_season(i,l) = nanmean(nanmean(double(temp_d) .* (area_mask./area_mask),1),2);
        case 'no3'
            clearvars temp_d
        temp_d = (squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) .* 1000) ; shading flat; 
        spatial_mean(i,l) = nanmean(nanmean(double(temp_d) .* (area_mask./area_mask),1),2);
        case 'o2'
            clearvars temp_d
        temp_d = (squeeze(mean(d_mean(:,:,(i-1)*3+1:i*3),3)) .* 1000) ; shading flat; 
        spatial_mean(i,l) = nanmean(nanmean(double(temp_d) .* (area_mask./area_mask),1),2);
    end
    
end
%  nanmax(nanmax(chl_mean .* 10^6))
%  nanmin(nanmin(chl_mean .* 10^6))

switch(var_name)
 case 'chl'
     clearvars temp_d
     for iij=1:12
        temp_d = (squeeze(d_mean(:,:,iij)) .* 10^6) ; shading flat; %kg/m^3 to ug/L
        spatial_mean(iij,l) = nanmean(nanmean(double(temp_d) .* (area_mask./area_mask),1),2);
     end
 end

end
      
save([var_name ,'_interp_historical_clim_spmean_', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'.mat']);

figure;
for l = 1:length(mod_name_pre)
hold on; plot(spatial_mean(:,l));
end
legend(mod_name_pre)

%% Ensemble
close all; clear; clc;

% area_mask=NaN(size(lon2,1),size(lon2,2));
% area_mask(find(lon2 >= 126.5 & lon2 <= 130 & 32.5<= lat2 & 35.4 >= lat2))=1;
% area_mask(index_nan) = NaN;
% extract_other_seas=find(isnan(area_mask)==1);
% monthly_mean(extract_other_seas) = NaN;

% time_dyear = 1965:2005; 
time_dyear = 1997:2005; 
var_name = 'chl'
% var_name = 'no3'
% var_name = 'o2'
mod_name_pre = {'IPSL-CM5A-LR','IPSL-CM5A-MR','IPSL-CM5B-LR', ...
    'GFDL-ESM2G','CESM1-BGC'}; 
% mod_name_pre = {'CNRM-CM5'}; 

% standard_name = 'mole_concentration_of_nitrate_in_sea_water'
% long_name     = 'Dissolved Nitrate Concentration at Surface'
% units         = 'mol m-3'
% 1 mol/m3 = 1000 mmol/m3

x=ncread([var_name '_merged_CESM1_BGC_his_1965-2005.nc'],'lon');
y=ncread([var_name '_merged_CESM1_BGC_his_1965-2005.nc'],'lat');

[lat lon]=meshgrid(y,x);
ik = 0;
for l = 1:length(mod_name_pre)
    clearvars -except time_dyear var_name mod_name_pre l lon lat meg_d ik
ik=ik+1;
mod_name = mod_name_pre{l};
disp(mod_name)
%     D:\장기생태\Dynamic\CMIP5_bio\no3
cd(['D:\장기생태\Dynamic\CMIP5_bio\' var_name])

temp_d=load([var_name ,'_interp_' mod_name '_historical_clim_', num2str(time_dyear(1)), '-', num2str(time_dyear(end)), '.mat'],'d_mean');
if ik == 1
    meg_d=[temp_d.d_mean];
else
    meg_d=(meg_d + temp_d.d_mean);
end
end
meg_data=meg_d ./ length(mod_name_pre);

seasonal = {'winter','spring','summer','fall'};
for i = 1:4   
%     d_mean(d_mean == 0) = NaN;

    fig=figure;
    switch(var_name)
        case 'chl'
        pcolor(lon,lat,squeeze(mean(meg_data(:,:,(i-1)*3+1:i*3),3)) .* 10^6); shading flat; 
        h=colorbar; caxis([0 10]);
        colorbarlabel='Chl (ug/L)'
        case 'no3'
        pcolor(lon,lat,squeeze(mean(meg_data(:,:,(i-1)*3+1:i*3),3)) .* 1000); shading flat; 
        h=colorbar; caxis([0 10]);
        colorbarlabel='no_3 (mmol/m^3)';
        case 'o2'
        pcolor(lon,lat,squeeze(mean(meg_data(:,:,(i-1)*3+1:i*3),3)) .* 1000); shading flat;      
        h=colorbar; caxis([0 380]);
        colorbarlabel='O_2 (mmol/m^3)';
    end
    title([var_name ,' interp ens historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' seasonal{i}]);
    ylim([15 52.0063])
    xlim([115 161.90])
    plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
    ylim([15 52.0063])
    xlim([115 161.90])
    colormap('jet'); set(get(h,'title'),'string',colorbarlabel);
    saveas(fig,[var_name ,' interp ens historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'. ' seasonal{i} '.png']);
    close;
end

sample=mean(meg_data,3);


fig=figure;
    switch(var_name)
        case 'chl'
        pcolor(lon,lat,sample .* 10^6); shading flat;  
        h=colorbar; caxis([0 10]);
        colorbarlabel='Chl (ug/L)';
        case 'no3'
        pcolor(lon,lat,sample .* 1000); shading flat;  
        h=colorbar; caxis([0 10]);
        colorbarlabel='no_3 (mmol/m^3)';
        case 'o2'
        pcolor(lon,lat,sample .* 1000); shading flat;      
        h=colorbar; caxis([0 380]);
        colorbarlabel='O_2 (mmol/m^3)';
    end 
    title([var_name ,' interp ens historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end))]);
    ylim([15 52.0063])
    xlim([115 161.90])
    plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
    ylim([15 52.0063])
    xlim([115 161.90])
%     colormap('jet')
    set(get(h,'title'),'string',colorbarlabel);
    saveas(fig,[var_name ,' interp ens historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'_yr.png']);
    close;
    
    
    fig=figure;
    pcolor(lon,lat,sample .* 1000); shading flat; 
    colorbar; 
    caxis([0 10])
    title([var_name ,' interp ens historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end))]);
    ylim([30 45])
    xlim([115 145])
    plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
    ylim([30 45])
    xlim([115 145])
%     colormap('jet')
    saveas(fig,[var_name ,' interp ens historical clim ', num2str(time_dyear(1)), '-', num2str(time_dyear(end)),'_yr.png']);
    close;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;
% file_SODA = '/home/cshwa/WOA/WOA05_NITR_00an1.nc';
% file_SODA_p = '/home/cshwa/WOA/WOA05_PHOS_00an1.nc';

% var_name = 'o2'
var_name = 'o2'

% n_name='woa05_no3_an_seasonal.nc';
n_name=['woa05_' var_name '_an_seasonal.nc'];

lon_tgt = double( ncread( n_name, 'X' ));
lat_tgt = double( ncread( n_name, 'Y' ));
% lev_tgt = double( ncread( file_SODA, 'depth' ));
switch(var_name)
    case 'o2'
    data_tgt = double( ncread( n_name, 'O2' )); %nitrate
    case 'no3'
    data_tgt = double( ncread( n_name, 'nitrate' )); %nitrate
end
% po4_tgt = double( ncread( file_SODA_p, 'p00an1' ));
% 
% n_lev_tgt = length(lev_tgt);
[lon, lat] = ndgrid( lon_tgt, lat_tgt );

data_p=squeeze(data_tgt(:,:,1,:)); % layer

data_p(data_p < -90) = NaN;

area_mask=NaN(size(lon,1),size(lon,2));
area_mask(find(lon >= 126.5 & lon <= 130 & 32.5<= lat & 35.4 >= lat))=1;
seasonal = {'winter','spring','summer','fall'};

for i = 1:4  %winter, spring, summer, fall (1~3, 4~6, 7~9, 10~12);
fig=figure;
switch(var_name)
        case 'no3'
        pcolor(lon,lat,squeeze(data_p(:,:,i))); shading flat;  %no3
        h=colorbar; caxis([0 10]);
        colorbarlabel='no_3 (mmol/m^3)';
        case 'o2'
        pcolor(lon,lat,squeeze(data_p(:,:,i))  .* 10^3./22.391); shading flat;  %oxygen ml/l to mmol/m3      
        h=colorbar; caxis([0 380]);
        colorbarlabel='O_2 (mmol/m^3)';
    end 

%     title(['no3 WOA05 clim monthly ', num2str(i), '-mon']);
    title([var_name ' WOA05 clim monthly ', seasonal{i}]);
    ylim([15 52.0063])
    xlim([115 161.90])
    plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
    ylim([15 52.0063])
    xlim([115 161.90])
    colormap('jet'); set(get(h,'title'),'string',colorbarlabel);
%     saveas(fig,['no3_woa05_clim_monthly_', num2str(i), '-mon.png']);
    saveas(fig,[var_name '_woa05_clim_monthly_', seasonal{i}, '.png']);
    close;
    
   %%  spatial mean
%     fig=figure;    
    switch(var_name)
        case 'no3'
        temp_d = squeeze(data_p(:,:,i));
        spatial_mean(i) = nanmean(nanmean(double(temp_d) .* (area_mask./area_mask),1),2);
        case 'o2'
        temp_d = (squeeze(data_p(:,:,i))  .* 10^3./22.391);
        spatial_mean(i) = nanmean(nanmean(double(temp_d) .* (area_mask./area_mask),1),2);
    end 
end
save([var_name ,'_WOA05_spmean.mat']);

%% obs vs. model time serise
close all; clear; clc;
woa_no3=load('D:\장기생태\Dynamic\CMIP5_bio\WOA\seasonal\no3_WOA05_spmean.mat');
woa_o2=load('D:\장기생태\Dynamic\CMIP5_bio\WOA\seasonal\o2_WOA05_spmean.mat');
sat_chl=load('D:\장기생태\Dynamic\CMIP5_bio\sat\seaWIFS_spmean_monthly.mat');
mod_no3 = load('D:\장기생태\Dynamic\CMIP5_bio\no3\no3_interp_historical_clim_spmean_1965-2005.mat');
mod_o2 = load('D:\장기생태\Dynamic\CMIP5_bio\o2\o2_interp_historical_clim_spmean_1965-2005.mat');
mod_chl = load('D:\장기생태\Dynamic\CMIP5_bio\chl\chl_interp_historical_clim_spmean_1997-2005.mat');

figure;
for l = 1:length(mod_o2.mod_name_pre)
hold on; plot(mod_o2.spatial_mean(:,l));
end
plot(woa_o2.spatial_mean,'r--','linew',2);
legend([mod_o2.mod_name_pre, 'WOA05']);
xticks(1:4); xticklabels({'winter','spring','summer','fall'}); grid on;
set(gca,'fontsize',13,'fontweight','bold');
ylabel('O2 (mmol/m^3)');

figure;
for l = 1:length(mod_no3.mod_name_pre)
hold on; plot(mod_no3.spatial_mean(:,l),'linew',2);
end
plot(woa_no3.spatial_mean,'r--','linew',2);
legend([mod_no3.mod_name_pre, 'WOA05']);
xticks(1:4); xticklabels({'winter','spring','summer','fall'}); grid on;
set(gca,'fontsize',13,'fontweight','bold');
ylabel('NO3 (mmol/m^3)');

  


