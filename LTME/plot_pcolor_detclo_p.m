close all; clear; clc; 

grd_file='grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');

mod_file='D:\장기생태\Dynamic\result\2001\2001_mp_p_detclo.nc';
mod_no3=ncread(mod_file, 'NO3');
mod_nh4=ncread(mod_file, 'NH4');
mod_chl=ncread(mod_file, 'chlorophyll');
mod_do=ncread(mod_file, 'oxygen');
mod_temp=ncread(mod_file, 'temp');
mod_salt=ncread(mod_file, 'salt');
mod_po4=ncread(mod_file, 'tPO4');
mod_LdetP=ncread(mod_file, 'LdetritusP');
mod_LdetN=ncread(mod_file, 'LdetritusN');
mod_SdetP=ncread(mod_file, 'SdetritusP');
mod_SdetN=ncread(mod_file, 'SdetritusN');

% make 2001 monthly mean
eom_d = 0
k=0
for i = 2001:2001
            k=k+1;
    for j = 1:12
        if j==1
            eom_d(k,j) = eom_d(k,1) + eomday(i,j); % 1980 is leap-yr
        else 
            eom_d(k,j) = eom_d(k,j-1) + eomday(i,j); % 1980 is leap-yr
        end
    end
end


for i = 1:12
if i==1
    LdetP_mon(:,:,:,i) = squeeze(nanmean(mod_LdetP(:,:,:,1 : eom_d(1,i)),4));
    SdetP_mon(:,:,:,i) = squeeze(nanmean(mod_SdetP(:,:,:,1 : eom_d(1,i)),4));
    LdetN_mon(:,:,:,i) = squeeze(nanmean(mod_LdetN(:,:,:,1 : eom_d(1,i)),4));
    SdetN_mon(:,:,:,i) = squeeze(nanmean(mod_SdetN(:,:,:,1 : eom_d(1,i)),4));
    no3_mon(:,:,:,i) = squeeze(nanmean(mod_no3(:,:,:,1 : eom_d(1,i)),4));
    nh4_mon(:,:,:,i) = squeeze(nanmean(mod_nh4(:,:,:,1 : eom_d(1,i)),4));
    po4_mon(:,:,:,i) = squeeze(nanmean(mod_po4(:,:,:,1 : eom_d(1,i)),4));
    chl_mon(:,:,:,i) = squeeze(nanmean(mod_chl(:,:,:,1 : eom_d(1,i)),4));
else
    LdetP_mon(:,:,:,i) = squeeze(nanmean(mod_LdetP(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
    SdetP_mon(:,:,:,i) = squeeze(nanmean(mod_SdetP(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
    LdetN_mon(:,:,:,i) = squeeze(nanmean(mod_LdetN(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
    SdetN_mon(:,:,:,i) = squeeze(nanmean(mod_SdetN(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
    no3_mon(:,:,:,i) = squeeze(nanmean(mod_no3(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
    nh4_mon(:,:,:,i) = squeeze(nanmean(mod_nh4(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
    po4_mon(:,:,:,i) = squeeze(nanmean(mod_po4(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
    chl_mon(:,:,:,i) = squeeze(nanmean(mod_chl(:,:,:,eom_d(1,i-1)+1 : eom_d(1,i)),4));
end
end

s_lv = 20; %surface
b_lv = 1; %bottom

%% Large det.

clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(LdetP_mon(:,:,s_lv,:))))):0.01:nanmax(nanmax(nanmax(squeeze(LdetP_mon(:,:,s_lv,:)))))
for i =1:12    
fig = figure; hold on;
pcolor(x,y,squeeze(LdetP_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL LdetP ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','LdetP(uM)','fontsize',14,'fontweight','bold');
caxis([0 0.1]);

[cs,h]=contour(x,y,squeeze(LdetP_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_LdetP_' num2str(i) '-mon'],'-dpng')
end

clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(LdetN_mon(:,:,s_lv,:))))):0.1:nanmax(nanmax(nanmax(squeeze(LdetN_mon(:,:,s_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(LdetN_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(s_lv) 'lv LdetN ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','LdetN(uM)','fontsize',14,'fontweight','bold');
caxis([0 1.8]);

[cs,h]=contour(x,y,squeeze(LdetN_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_LdetN_' num2str(i) '-mon'],'-dpng')
end

%bot
clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(LdetP_mon(:,:,b_lv,:))))):0.01:nanmax(nanmax(nanmax(squeeze(LdetP_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(LdetP_mon(:,:,b_lv,i))); shading flat;  
contour(x,y,squeeze(LdetP_mon(:,:,b_lv,i)),'k')
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL LdetP ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','LdetP(uM)','fontsize',14,'fontweight','bold');
caxis([0 0.1]);

[cs,h]=contour(x,y,squeeze(LdetP_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_LdetP_' num2str(i) '-mon'],'-dpng')
end

clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(LdetN_mon(:,:,b_lv,:))))):0.1:nanmax(nanmax(nanmax(squeeze(LdetN_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(LdetN_mon(:,:,b_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(b_lv) 'lv LdetN ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','LdetN(uM)','fontsize',14,'fontweight','bold');
% caxis([0 0.1]);
caxis([0 3]);
[cs,h]=contour(x,y,squeeze(LdetN_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_LdetN_' num2str(i) '-mon'],'-dpng')
end

%% Small det.

clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(SdetP_mon(:,:,s_lv,:))))):0.01:nanmax(nanmax(nanmax(squeeze(SdetP_mon(:,:,s_lv,:)))))
for i =1:12    
fig = figure; hold on;
pcolor(x,y,squeeze(SdetP_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL SdetP ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','SdetP(uM)','fontsize',14,'fontweight','bold');
caxis([0 0.3]);

[cs,h]=contour(x,y,squeeze(SdetP_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_SdetP_' num2str(i) '-mon'],'-dpng')
end

clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(SdetN_mon(:,:,s_lv,:))))):0.2:nanmax(nanmax(nanmax(squeeze(SdetN_mon(:,:,s_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(SdetN_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(s_lv) 'lv SdetN ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','SdetN(uM)','fontsize',14,'fontweight','bold');
caxis([0 6]);

[cs,h]=contour(x,y,squeeze(SdetN_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_SdetN_' num2str(i) '-mon'],'-dpng')
end

%bot
clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(SdetP_mon(:,:,b_lv,:))))):0.01:nanmax(nanmax(nanmax(squeeze(SdetP_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(SdetP_mon(:,:,b_lv,i))); shading flat;  
contour(x,y,squeeze(SdetP_mon(:,:,b_lv,i)),'k')
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL SdetP ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','SdetP(uM)','fontsize',14,'fontweight','bold');
caxis([0 0.2]);

[cs,h]=contour(x,y,squeeze(SdetP_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_SdetP_' num2str(i) '-mon'],'-dpng')
end

clearvars temp_range
temp_range=nanmin(nanmin(nanmin(squeeze(SdetN_mon(:,:,b_lv,:))))):0.2:nanmax(nanmax(nanmax(squeeze(SdetN_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(SdetN_mon(:,:,b_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(b_lv) 'lv SdetN ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','SdetN(uM)','fontsize',14,'fontweight','bold');
caxis([0 6]);
[cs,h]=contour(x,y,squeeze(SdetN_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_SdetN_' num2str(i) '-mon'],'-dpng')
end

%% nut.

clearvars temp_range
temp_range=0:0.5:nanmax(nanmax(nanmax(squeeze(no3_mon(:,:,s_lv,:)))))
for i =1:12    
fig = figure; hold on;
pcolor(x,y,squeeze(no3_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL no3 ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','no3(uM)','fontsize',14,'fontweight','bold');
caxis([0 15]);

[cs,h]=contour(x,y,squeeze(no3_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_no3_' num2str(i) '-mon'],'-dpng')
end

clearvars temp_range
temp_range=0:0.5:nanmax(nanmax(nanmax(squeeze(nh4_mon(:,:,s_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(nh4_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(s_lv) 'lv nh4 ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','nh4(uM)','fontsize',14,'fontweight','bold');
caxis([0 10]);

[cs,h]=contour(x,y,squeeze(nh4_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_nh4_' num2str(i) '-mon'],'-dpng')
end

%bot
clearvars temp_range
temp_range=0:1:nanmax(nanmax(nanmax(squeeze(no3_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(no3_mon(:,:,b_lv,i))); shading flat;  
contour(x,y,squeeze(no3_mon(:,:,b_lv,i)),'k')
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL no3 ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','no3(uM)','fontsize',14,'fontweight','bold');
caxis([0 15]);

[cs,h]=contour(x,y,squeeze(no3_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_no3_' num2str(i) '-mon'],'-dpng')
end

clearvars temp_range
temp_range=0:0.5:nanmax(nanmax(nanmax(squeeze(nh4_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(nh4_mon(:,:,b_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(b_lv) 'lv nh4 ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','nh4(uM)','fontsize',14,'fontweight','bold');
% caxis([0 0.1]);
caxis([0 10]);
[cs,h]=contour(x,y,squeeze(nh4_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_nh4_' num2str(i) '-mon'],'-dpng')
end

%% PO4 
clearvars temp_range
temp_range=0:0.2:nanmax(nanmax(nanmax(squeeze(po4_mon(:,:,s_lv,:)))))
for i =1:12    
fig = figure; hold on;
pcolor(x,y,squeeze(po4_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL po4 ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','po4(uM)','fontsize',14,'fontweight','bold');
caxis([0 1]);

[cs,h]=contour(x,y,squeeze(po4_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_po4_' num2str(i) '-mon'],'-dpng')
end

%bot
clearvars temp_range
temp_range=0:0.1:nanmax(nanmax(nanmax(squeeze(po4_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(po4_mon(:,:,b_lv,i))); shading flat;  
contour(x,y,squeeze(po4_mon(:,:,b_lv,i)),'k')
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL po4 ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','po4(uM)','fontsize',14,'fontweight','bold');
caxis([0 1]);

[cs,h]=contour(x,y,squeeze(po4_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_po4_' num2str(i) '-mon'],'-dpng')
end

%% chl
clearvars temp_range
temp_range=0:0.5:nanmax(nanmax(nanmax(squeeze(chl_mon(:,:,s_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(chl_mon(:,:,s_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(s_lv) 'lv chl ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','chl(uM)','fontsize',14,'fontweight','bold');
% caxis([0 0.1]);
caxis([0 10]);
[cs,h]=contour(x,y,squeeze(chl_mon(:,:,s_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(s_lv) 'lv_MODEL_chl_' num2str(i) '-mon'],'-dpng')
end


clearvars temp_range
temp_range=0:0.5:nanmax(nanmax(nanmax(squeeze(chl_mon(:,:,b_lv,:)))))
for i =1:12
fig = figure; hold on;
pcolor(x,y,squeeze(chl_mon(:,:,b_lv,i))); shading flat;        
ylim([min(min(y)) max(max(y))])
xlim([min(min(x)) max(max(x))])
title(['GY 2001 MODEL ' num2str(b_lv) 'lv chl ' num2str(i) '-mon']);
xlabel('lon','fontsize',16)
ylabel('lat','fontsize',16)
h = colorbar; set(get(h,'title'),'string','chl(uM)','fontsize',14,'fontweight','bold');
% caxis([0 0.1]);
caxis([0 10]);
[cs,h]=contour(x,y,squeeze(chl_mon(:,:,b_lv,i)),temp_range,'linecolor',[.01 .01 .01],'linewidth',1);
clabel(cs,h,'fontsize',14,'labelspacing',700,'fontweight','bold');
grid on;
set(gca,'fontsize',13,'fontweight','bold')
print(fig,['2001_' num2str(b_lv) 'lv_MODEL_chl_' num2str(i) '-mon'],'-dpng')
end


