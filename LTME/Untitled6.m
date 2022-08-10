clearvars salt_clim_*
k=0;
for i = [2,5,8,11]   
    k=k+1;
            salt_clim_surf(:,:,k) = squeeze(mean(salt_int_surf(:,:,i:12:end),3));
            salt_clim_bot(:,:,k) = squeeze(mean(salt_int_bot(:,:,i:12:end),3));
            salt_clim_surf_1st(:,:,k) = squeeze(mean(salt_int_surf_1st(:,:,i:12:end),3));
            salt_clim_bot_1st(:,:,k) = squeeze(mean(salt_int_bot_1st(:,:,i:12:end),3));
            salt_clim_surf_2nd(:,:,k) = squeeze(mean(salt_int_surf_2nd(:,:,i:12:end),3));
            salt_clim_bot_2nd(:,:,k) = squeeze(mean(salt_int_bot_2nd(:,:,i:12:end),3));           
end

fig=figure; 
pcolor(x,y,squeeze(salt_clim_surf(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title('1997-2019 surface salinity OBS','fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])
caxis([30 34])
saveas(fig,['all_obs_surf_salt.png']);

fig=figure; 
pcolor(x,y,squeeze(salt_clim_surf_1st(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title('2007-2015 surface salinity OBS','fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])
caxis([30 34])
saveas(fig,['1st_obs_surf_salt.png']);

fig=figure; 
pcolor(x,y,squeeze(salt_clim_surf_2nd(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title('2016-2019 surface salinity OBS','fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])
caxis([30 34])
saveas(fig,['2nd_obs_surf_salt.png']);

fig=figure; 
pcolor(x,y,squeeze(salt_clim_bot(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_bot_08(:,:,2)).*(mask./mask));
title('1997-2019 bottom salinity OBS','fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])
caxis([30 34])
saveas(fig,['all_obs_bot_salt.png']);

fig=figure; 
pcolor(x,y,squeeze(salt_clim_bot_1st(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_bot_08(:,:,2)).*(mask./mask));
title('2007-2015 bottom salinity OBS','fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])
caxis([30 34])
saveas(fig,['1st_obs_bot_salt.png']);

fig=figure; 
pcolor(x,y,squeeze(salt_clim_bot_2nd(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_bot_08(:,:,2)).*(mask./mask));
title('2016-2019 bottom salinity OBS','fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])
caxis([30 34])
saveas(fig,['2nd_obs_bot_salt.png']);
