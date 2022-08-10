close all; clear; clc;

yoonja=load('D:\Àå±â»ýÅÂ\Dynamic\06_river\yoonjakangs_koem_data_monthly_v2_16points.mat'); % yoonjakangs_koem_data_processing_v2_16points.m
[res I]=sort([1,5,4,3,2,6,8,9,7]);
% sorting confirm
nn=['±¤¾çÇ×'; '±¤¾ç4'; '±¤¾ç3'; '±¤¾ç2'; '±¤¾ç1'; '±¤¾ç5'; '¿©¼ö2'; '¿©¼ö3'; '¿©¼ö1'];

% combining the tag and outter point excluding
name_tag = nn(I,:); 
size_tag = length(name_tag);

% sample st. num
sp_9p_st = [1:6, 14:16]; 
sp_8p_st = [2:6, 14:16];
sp_12p_st = [2:10, 14:16];
sp_13p_st = [1:10, 14:16];
sp_15p_st = [2:16];
sp_16p_st = [1:16];

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for l = 2:13
    k=k+1;
name_tag{l}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for l = 14:16
    k=k+1;
name_tag{l}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end




% location
lat_st = [34.8625,34.85139,34.88389,34.90222,34.92083,34.83194,34.86944,...
    34.90278,34.89167,34.9025,34.86361,34.91,34.94,34.73556,34.76278,34.76444];
lon_st = [127.7403,127.6797,127.6506,127.6822,127.8233,127.8011,127.7889,...
    127.7989,127.7597,127.7239,127.7103,127.7,127.8264,127.7661,127.7614,...
    127.8053];

refyear=2021;

clearvars lat_1 lon_1 lat_2 lat_3 lon_2 lon_3
grd_file='D:\Àå±â»ýÅÂ\Dynamic\KOEM\grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
h=ncread(grd_file, 'h');


temp_range=[-10:0.2:10]; temp_range_c=[-1:0.1:1];
salt_range=[-10:0.5:10]; salt_range_c=[-3:0.2:3];
chl_range=[-10:.5:10]; chl_range_c=[-3:.2:3];
no3_range=[-30:2:30]; no3_range_c=[-10:.5:10];
nh4_range=[-10:0.1:10]; nh4_range_c=[-1:.1:1]; 
po4_range=[0:.1:17]; po4_range_c=[0:.2:3];
phy_range=[0:.1:5]; phy_range_c=[0:1:3.5];
zoo_range=[0:.1:6]; zoo_range_c=[0:.1:1];
z_range=[0:.01:1]; z_range_c=[0:.01:0.3];


%% too coarse to provide fine pcolor so remove over heated temp 
% 
% % loc_ff=[[lon_obs(13)-0.001, lat_obs(13)+0.001];[lon_obs(13)-0.001, lat_obs(13)-0.001]; [lon_obs(13)+0.001, lat_obs(13)+0.001]; [lon_obs(13)+0.001, lat_obs(13)-0.001]];
% % loc_ff=[[lon_obs(13)-0.01, lat_obs(13)+0.01];[lon_obs(13)-0.01, lat_obs(13)-0.01]; [lon_obs(13)+0.01, lat_obs(13)+0.01]; [lon_obs(13)+0.01, lat_obs(13)-0.01]];
% loc_ff=[[lon_obs(13)-0.005, lat_obs(13)+0.005];[lon_obs(13)-0.005, lat_obs(13)-0.005]; [lon_obs(13)+0.005, lat_obs(13)+0.005]; [lon_obs(13)+0.005, lat_obs(13)-0.005]];
% 
% lon_obs_f = [lon_obs; loc_ff(:,1)];
% lat_obs_f = [lat_obs; loc_ff(:,2)];
% temp_obs_surf_f = [temp_obs_surf; repmat(temp_obs_surf(4),4,1)];

%% 2008, 2013, 2015, 2016

salt_sur_08=yoonja.salt_sur(:,length(1997:2008 - 1)*12+1:length(1997:2008)*12);
salt_bot_08=yoonja.salt_bot(:,length(1997:2008 - 1)*12+1:length(1997:2008)*12);

salt_sur_13=yoonja.salt_sur(:,length(1997:2013 - 1)*12+1:length(1997:2013)*12);
salt_bot_13=yoonja.salt_bot(:,length(1997:2013 - 1)*12+1:length(1997:2013)*12);

salt_sur_15=yoonja.salt_sur(:,length(1997:2015 - 1)*12+1:length(1997:2015)*12);
salt_bot_15=yoonja.salt_bot(:,length(1997:2015 - 1)*12+1:length(1997:2015)*12);

salt_sur_16=yoonja.salt_sur(:,length(1997:2016 - 1)*12+1:length(1997:2016)*12);
salt_bot_16=yoonja.salt_bot(:,length(1997:2016 - 1)*12+1:length(1997:2016)*12);

salt_sur_17=yoonja.salt_sur(:,length(1997:2017 - 1)*12+1:length(1997:2017)*12);
salt_bot_17=yoonja.salt_bot(:,length(1997:2017 - 1)*12+1:length(1997:2017)*12);

salt_sur_18=yoonja.salt_sur(:,length(1997:2018 - 1)*12+1:length(1997:2018)*12);
salt_bot_18=yoonja.salt_bot(:,length(1997:2018 - 1)*12+1:length(1997:2018)*12);

salt_sur_19=yoonja.salt_sur(:,length(1997:2019 - 1)*12+1:length(1997:2019)*12);
salt_bot_19=yoonja.salt_bot(:,length(1997:2019 - 1)*12+1:length(1997:2019)*12);


% linear interpolation
% temp_int_surf=griddata(lon_obs, lat_obs,temp_obs_surf,lon,lat);

clearvars pre_*
for i = 1:12
clearvars pre_*
if i == 5 || i == 11
    disp(i)
pre_salt_int_surf_08=griddata(lon_st(sp_8p_st), lat_st(sp_8p_st),squeeze(salt_sur_08((sp_8p_st),i)),x,y);
salt_int_surf_08(:,:,i)=pre_salt_int_surf_08;
pre_salt_int_bot_08=griddata(lon_st(sp_8p_st), lat_st(sp_8p_st),squeeze(salt_bot_08((sp_8p_st),i)),x,y);
salt_int_bot_08(:,:,i)=pre_salt_int_bot_08;

pre_salt_int_surf_13=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_sur_13((sp_15p_st),i)),x,y);
salt_int_surf_13(:,:,i)=pre_salt_int_surf_13;
pre_salt_int_bot_13=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_bot_13((sp_15p_st),i)),x,y);
salt_int_bot_13(:,:,i)=pre_salt_int_bot_13;

pre_salt_int_surf_15=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_sur_15((sp_15p_st),i)),x,y);
salt_int_surf_15(:,:,i)=pre_salt_int_surf_15;
pre_salt_int_bot_15=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_bot_15((sp_15p_st),i)),x,y);
salt_int_bot_15(:,:,i)=pre_salt_int_bot_15;

pre_salt_int_surf_16=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_sur_16((sp_15p_st),i)),x,y);
salt_int_surf_16(:,:,i)=pre_salt_int_surf_16;
pre_salt_int_bot_16=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_bot_16((sp_15p_st),i)),x,y);
salt_int_bot_16(:,:,i)=pre_salt_int_bot_16;

pre_salt_int_surf_17=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_sur_17((sp_15p_st),i)),x,y);
salt_int_surf_17(:,:,i)=pre_salt_int_surf_17;
pre_salt_int_bot_17=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_bot_17((sp_15p_st),i)),x,y);
salt_int_bot_17(:,:,i)=pre_salt_int_bot_17;

pre_salt_int_surf_18=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_sur_18((sp_15p_st),i)),x,y);
salt_int_surf_18(:,:,i)=pre_salt_int_surf_18;
pre_salt_int_bot_18=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_bot_18((sp_15p_st),i)),x,y);
salt_int_bot_18(:,:,i)=pre_salt_int_bot_18;

pre_salt_int_surf_19=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_sur_19((sp_15p_st),i)),x,y);
salt_int_surf_19(:,:,i)=pre_salt_int_surf_19;
pre_salt_int_bot_19=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(salt_bot_19((sp_15p_st),i)),x,y);
salt_int_bot_19(:,:,i)=pre_salt_int_bot_19;

else
pre_salt_int_surf_08=griddata(lon_st(sp_9p_st), lat_st(sp_9p_st),squeeze(salt_sur_08((sp_9p_st),i)),x,y);
salt_int_surf_08(:,:,i)=pre_salt_int_surf_08;
pre_salt_int_bot_08=griddata(lon_st(sp_9p_st), lat_st(sp_9p_st),squeeze(salt_bot_08((sp_9p_st),i)),x,y);
salt_int_bot_08(:,:,i)=pre_salt_int_bot_08;

pre_salt_int_surf_13=griddata(lon_st, lat_st,squeeze(salt_sur_13(:,i)),x,y);
salt_int_surf_13(:,:,i)=pre_salt_int_surf_13;
pre_salt_int_bot_13=griddata(lon_st, lat_st,squeeze(salt_bot_13(:,i)),x,y);
salt_int_bot_13(:,:,i)=pre_salt_int_bot_13;

pre_salt_int_surf_15=griddata(lon_st, lat_st,squeeze(salt_sur_15(:,i)),x,y);
salt_int_surf_15(:,:,i)=pre_salt_int_surf_15;
pre_salt_int_bot_15=griddata(lon_st, lat_st,squeeze(salt_bot_15(:,i)),x,y);
salt_int_bot_15(:,:,i)=pre_salt_int_bot_15;

pre_salt_int_surf_16=griddata(lon_st, lat_st,squeeze(salt_sur_16(:,i)),x,y);
salt_int_surf_16(:,:,i)=pre_salt_int_surf_16;
pre_salt_int_bot_16=griddata(lon_st, lat_st,squeeze(salt_bot_16(:,i)),x,y);
salt_int_bot_16(:,:,i)=pre_salt_int_bot_16;

pre_salt_int_surf_17=griddata(lon_st, lat_st,squeeze(salt_sur_17(:,i)),x,y);
salt_int_surf_17(:,:,i)=pre_salt_int_surf_17;
pre_salt_int_bot_17=griddata(lon_st, lat_st,squeeze(salt_bot_17(:,i)),x,y);
salt_int_bot_17(:,:,i)=pre_salt_int_bot_17;

pre_salt_int_surf_18=griddata(lon_st, lat_st,squeeze(salt_sur_18(:,i)),x,y);
salt_int_surf_18(:,:,i)=pre_salt_int_surf_18;
pre_salt_int_bot_18=griddata(lon_st, lat_st,squeeze(salt_bot_18(:,i)),x,y);
salt_int_bot_18(:,:,i)=pre_salt_int_bot_18;

pre_salt_int_surf_19=griddata(lon_st, lat_st,squeeze(salt_sur_19(:,i)),x,y);
salt_int_surf_19(:,:,i)=pre_salt_int_surf_19;
pre_salt_int_bot_19=griddata(lon_st, lat_st,squeeze(salt_bot_19(:,i)),x,y);
salt_int_bot_19(:,:,i)=pre_salt_int_bot_19;
end

% pre_salt_int_surf_08=kriging(lon_st(sp_9p_st)', lat_st(sp_9p_st)',squeeze(salt_sur_08((sp_9p_st),i)),x,y); % x, y would be row vector
% salt_int_surf_08(:,:,i)=pre_salt_int_surf_08;
% pre_salt_int_bot_08=kriging(lon_st(sp_9p_st)', lat_st(sp_9p_st)',squeeze(salt_bot_08((sp_9p_st),i)),x,y);
% salt_int_bot_08(:,:,i)=pre_salt_int_bot_08;
% 
% pre_salt_int_surf_13=kriging(lon_st', lat_st',squeeze(salt_sur_13(:,i)),x,y);
% salt_int_surf_13(:,:,i)=pre_salt_int_surf_13;
% pre_salt_int_bot_13=kriging(lon_st', lat_st',squeeze(salt_bot_13(:,i)),x,y);
% salt_int_bot_13(:,:,i)=pre_salt_int_bot_13;
% 
% pre_salt_int_surf_15=kriging(lon_st', lat_st',squeeze(salt_sur_15(:,i)),x,y);
% salt_int_surf_15(:,:,i)=pre_salt_int_surf_15;
% pre_salt_int_bot_15=kriging(lon_st', lat_st',squeeze(salt_bot_15(:,i)),x,y);
% salt_int_bot_15(:,:,i)=pre_salt_int_bot_15;
% 
% pre_salt_int_surf_16=kriging(lon_st', lat_st',squeeze(salt_sur_16(:,i)),x,y);
% salt_int_surf_16(:,:,i)=pre_salt_int_surf_16;
% pre_salt_int_bot_16=kriging(lon_st', lat_st',squeeze(salt_bot_16(:,i)),x,y);
% salt_int_bot_16(:,:,i)=pre_salt_int_bot_16;
% 
% pre_salt_int_surf_17=kriging(lon_st', lat_st',squeeze(salt_sur_17(:,i)),x,y);
% salt_int_surf_17(:,:,i)=pre_salt_int_surf_17;
% pre_salt_int_bot_17=kriging(lon_st', lat_st',squeeze(salt_bot_17(:,i)),x,y);
% salt_int_bot_17(:,:,i)=pre_salt_int_bot_17;
% 
% pre_salt_int_surf_18=kriging(lon_st', lat_st',squeeze(salt_sur_18(:,i)),x,y);
% salt_int_surf_18(:,:,i)=pre_salt_int_surf_18;
% pre_salt_int_bot_18=kriging(lon_st', lat_st',squeeze(salt_bot_18(:,i)),x,y);
% salt_int_bot_18(:,:,i)=pre_salt_int_bot_18;
% 
% pre_salt_int_surf_19=kriging(lon_st', lat_st',squeeze(salt_sur_19(:,i)),x,y);
% salt_int_surf_19(:,:,i)=pre_salt_int_surf_19;
% pre_salt_int_bot_19=kriging(lon_st', lat_st',squeeze(salt_bot_19(:,i)),x,y);
% salt_int_bot_19(:,:,i)=pre_salt_int_bot_19;
end

clearvars pre_salt_int_surf salt_int_bot
for i = 1:length(yoonja.salt_sur)
    clearvars pre_salt_int_surf pre_salt_int_bot
    if i < 146
        if mod(i,12) == 5 || mod(i,12) == 11
            pre_salt_int_surf=griddata(lon_st(sp_8p_st), lat_st(sp_8p_st),squeeze(yoonja.salt_sur((sp_8p_st),i)),x,y);
            salt_int_surf(:,:,i)=pre_salt_int_surf;
            pre_salt_int_bot=griddata(lon_st(sp_8p_st), lat_st(sp_8p_st),squeeze(yoonja.salt_bot((sp_8p_st),i)),x,y);
            salt_int_bot(:,:,i)=pre_salt_int_bot;
        else
            pre_salt_int_surf=griddata(lon_st(sp_9p_st), lat_st(sp_9p_st),squeeze(yoonja.salt_sur((sp_9p_st),i)),x,y);
            salt_int_surf(:,:,i)=pre_salt_int_surf;
            pre_salt_int_bot=griddata(lon_st(sp_9p_st), lat_st(sp_9p_st),squeeze(yoonja.salt_bot((sp_9p_st),i)),x,y);
            salt_int_bot(:,:,i)=pre_salt_int_bot;
        end
    elseif i >= 146 && i < 194
        if mod(i,12) == 5 || mod(i,12) == 11
            pre_salt_int_surf=griddata(lon_st(sp_12p_st), lat_st(sp_12p_st),squeeze(yoonja.salt_sur((sp_12p_st),i)),x,y);
            salt_int_surf(:,:,i)=pre_salt_int_surf;
            pre_salt_int_bot=griddata(lon_st(sp_12p_st), lat_st(sp_12p_st),squeeze(yoonja.salt_bot((sp_12p_st),i)),x,y);
            salt_int_bot(:,:,i)=pre_salt_int_bot;
        else
            pre_salt_int_surf=griddata(lon_st(sp_13p_st), lat_st(sp_13p_st),squeeze(yoonja.salt_sur((sp_13p_st),i)),x,y);
            salt_int_surf(:,:,i)=pre_salt_int_surf;
            pre_salt_int_bot=griddata(lon_st(sp_13p_st), lat_st(sp_13p_st),squeeze(yoonja.salt_bot((sp_13p_st),i)),x,y);
            salt_int_bot(:,:,i)=pre_salt_int_bot;
        end 
    elseif i >= 194
        disp(i)
            if mod(i,12) == 5 || mod(i,12) == 11
            pre_salt_int_surf=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(yoonja.salt_sur((sp_15p_st),i)),x,y);
            salt_int_surf(:,:,i)=pre_salt_int_surf;
            pre_salt_int_bot=griddata(lon_st(sp_15p_st), lat_st(sp_15p_st),squeeze(yoonja.salt_bot((sp_15p_st),i)),x,y);
            salt_int_bot(:,:,i)=pre_salt_int_bot;
        else
            pre_salt_int_surf=griddata(lon_st(sp_16p_st), lat_st(sp_16p_st),squeeze(yoonja.salt_sur((sp_16p_st),i)),x,y);
            salt_int_surf(:,:,i)=pre_salt_int_surf;
            pre_salt_int_bot=griddata(lon_st(sp_16p_st), lat_st(sp_16p_st),squeeze(yoonja.salt_bot((sp_16p_st),i)),x,y);
            salt_int_bot(:,:,i)=pre_salt_int_bot;
            end
    end
end

% clearvars pre_salt_int_surf salt_int_bot
% for i = [2,5,8,11]
%     clearvars pre_salt_int_surf pre_salt_int_bot       
%             salt_clim_surf(:,:,i) = nanmean(salt_int_surf(:,:,i:12:end));
%             salt_clim_bot(:,:,i) = nanmean(salt_int_bot(:,:,i:12:end));
% end

salt_int_surf_0th = salt_int_surf(:,:,1:length(1997:2007 - 1)*12);
salt_int_bot_0th = salt_int_bot(:,:,1:length(1997:2007 - 1)*12);
salt_int_surf_1st = salt_int_surf(:,:,length(1997:2007 - 1)*12+1:length(1997:2015)*12);
salt_int_bot_1st = salt_int_bot(:,:,length(1997:2007 - 1)*12+1:length(1997:2015)*12);
salt_int_surf_2nd = salt_int_surf(:,:,length(1997:2016 - 1)*12+1:end);
salt_int_bot_2nd = salt_int_bot(:,:,length(1997:2016 - 1)*12+1:end);

return

% clearvars salt_clim_*
% k=0;
% for i = [2,5,8,11]   
%     k=k+1;
%             salt_clim_surf(:,:,k) = squeeze(nanmean(salt_int_surf(:,:,i:12:end),3));
%             salt_clim_bot(:,:,k) = squeeze(nanmean(salt_int_bot(:,:,i:12:end),3));
%             salt_clim_surf_0th(:,:,k) = squeeze(nanmean(salt_int_surf(:,:,i:12:end),3));
%             salt_clim_bot_0th(:,:,k) = squeeze(nanmean(salt_int_bot(:,:,i:12:end),3));
%             salt_clim_surf_1st(:,:,k) = squeeze(nanmean(salt_int_surf_1st(:,:,i:12:end),3));
%             salt_clim_bot_1st(:,:,k) = squeeze(nanmean(salt_int_bot_1st(:,:,i:12:end),3));
%             salt_clim_surf_2nd(:,:,k) = squeeze(nanmean(salt_int_surf_2nd(:,:,i:12:end),3));
%             salt_clim_bot_2nd(:,:,k) = squeeze(nanmean(salt_int_bot_2nd(:,:,i:12:end),3));           
% end

clearvars salt_clim_*
k=0;
for i = [2,5,8,11]   
    k=k+1;
            salt_clim_surf(:,:,k) = squeeze(mean(salt_int_surf(:,:,i:12:end),3));
            salt_clim_bot(:,:,k) = squeeze(mean(salt_int_bot(:,:,i:12:end),3));
            salt_clim_surf_0th(:,:,k) = squeeze(nanmean(salt_int_surf_0th(:,:,i:12:end),3));
            salt_clim_bot_0th(:,:,k) = squeeze(nanmean(salt_int_bot_0th(:,:,i:12:end),3));
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
pcolor(x,y,squeeze(salt_clim_surf_0th(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title('1997-2006 surface salinity OBS','fontsize',30);set(gca,'FontSize',16);
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
saveas(fig,['0th_obs_surf_salt.png']);

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
pcolor(x,y,squeeze(salt_clim_bot_0th(:,:,2)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title('1997-2006 bottom salinity OBS','fontsize',30);set(gca,'FontSize',16);
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
saveas(fig,['0th_obs_bot_salt.png']);

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



%% full time series
for i = 1997 : 2019
    for j = [2,5,8,11]
fig=figure; 
pcolor(x,y,squeeze(salt_int_surf(:,:,length(1997:i - 1)*12+j)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title([num2str(i),'-', num2str(j,'%02d') ,' surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars sp_st_confirm
if length(1997:i - 1)*12+j < 146
    sp_st_confirm = sp_9p_st;
elseif length(1997:i - 1)*12+j >= 146 && length(1997:i - 1)*12+j < 194
    sp_st_confirm = sp_13p_st;
elseif length(1997:i - 1)*12+j >= 194
    sp_st_confirm = sp_16p_st;
end

if mod(length(1997:i - 1)*12+j,12) == 5 || mod(length(1997:i - 1)*12+j,12) == 11
    sp_pick_st = sp_st_confirm(2:end);
else
    sp_pick_st = sp_st_confirm(1:end);
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 

for l=sp_pick_st
    if l== 11 | l== 15
        text(lon_st(l)-0.001,lat_st(l)-0.005, char(name_tag{l}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(l),lat_st(l),'.','color','k'); %plot only 97
    else
        text(lon_st(l)+0.003,lat_st(l), char(name_tag{l}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(l),lat_st(l),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,[num2str(i),'_' num2str(j,'%02d') '_obs_surf_salt.png']);

fig=figure; 
pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:i - 1)*12+j)).*(mask./mask))
% pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title([num2str(i),'-', num2str(j,'%02d') ,' bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars sp_st_confirm
if length(1997:i - 1)*12+j < 146
    sp_st_confirm = sp_9p_st;
elseif length(1997:i - 1)*12+j >= 146 && length(1997:i - 1)*12+j < 194
    sp_st_confirm = sp_13p_st;
elseif length(1997:i - 1)*12+j >= 194
    sp_st_confirm = sp_16p_st;
end

if mod(length(1997:i - 1)*12+j,12) == 5 || mod(length(1997:i - 1)*12+j,12) == 11
    sp_pick_st = sp_st_confirm(2:end);
else
    sp_pick_st = sp_st_confirm(1:end);
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 

for l=sp_pick_st
    if l== 11 | l== 15
        text(lon_st(l)-0.001,lat_st(l)-0.005, char(name_tag{l}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(l),lat_st(l),'.','color','k'); %plot only 97
    else
        text(lon_st(l)+0.003,lat_st(l), char(name_tag{l}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(l),lat_st(l),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,[num2str(i),'_' num2str(j,'%02d') '_obs_bot_salt.png']);

    end
end

%% 2008
cd D:\Àå±â»ýÅÂ\Dynamic\KOEM\koem_yoonja_kang\pic_hori

fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_surf_08(:,:,2)).*(mask./mask));
title(['2008 surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=[1:6,14:16]
    if i== 11 | i== 15
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2008_obs_surf_salt.png');


fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_bot_08(:,:,2)).*(mask./mask));
title(['2008 bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
% contourf(x,y,squeeze(salt_int_bot_08(:,:,2)).*(mask./mask))
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=[1:6,14:16]
    if i== 11 | i== 15
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
% caxis([30 33.8])
colormap jet
saveas(fig,'2008_obs_bot_salt.png');

%% else

% 2013
fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_surf_13(:,:,2)).*(mask./mask));
title(['2013 surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 15
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2013_obs_surf_salt.png');

fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_bot_13(:,:,2)).*(mask./mask));
title(['2013 bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 15
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2013_obs_bot_salt.png');

% 2015
fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_surf_15(:,:,2)).*(mask./mask));
title(['2015 surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 15
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2015_obs_surf_salt.png');

fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_bot_15(:,:,2)).*(mask./mask));
title(['2015 bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 15
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2015_obs_bot_salt.png');

% 2016
fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_surf_16(:,:,2)).*(mask./mask));
title(['2016 surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35]) 
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2016_obs_surf_salt.png');

fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_bot_16(:,:,2)).*(mask./mask));
title(['2016 bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2016_obs_bot_salt.png');

% 2017
fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_surf_17(:,:,2)).*(mask./mask));
title(['2017 surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2017_obs_surf_salt.png');

fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_bot_17(:,:,2)).*(mask./mask));
title(['2017 bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2017_obs_bot_salt.png');


% 2018
fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_surf_18(:,:,2)).*(mask./mask));
title(['2018 surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2018_obs_surf_salt.png');

fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_bot_18(:,:,2)).*(mask./mask));
title(['2018 bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2018_obs_bot_salt.png');

% 2019
fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_surf_19(:,:,2)).*(mask./mask));
title(['2019 surface salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2019_obs_surf_salt.png');

fig=figure; 
% pcolor(x,y,squeeze(salt_int_bot(:,:,length(1997:2013 - 1)*12+11)).*(mask./mask))
pcolor(x,y,squeeze(salt_int_bot_19(:,:,2)).*(mask./mask));
title(['2019 bottom salinity OBS'],'fontsize',30);set(gca,'FontSize',16);
shading flat
colorbar
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

clearvars name_tag
name_tag{1}={'±¤¾çÇ×01'};
k=0
for i = 2:13
    k=k+1;
name_tag{i}={['±¤¾ç¸¸',num2str(k,'%02d')]};
end
k=0
for i = 14:16
    k=k+1;
name_tag{i}={['¿©¼ö¿¬¾È',num2str(k,'%02d')]};
end

%plot only in the domain
% figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
% set(gca,'fontsize',16,'fontweight','bold'); 
grid on; 
for i=1:length(name_tag)
    if i== 11 | i== 16
        text(lon_st(i)-0.001,lat_st(i)-0.005, char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    else
        text(lon_st(i)+0.003,lat_st(i), char(name_tag{i}),'color','k','fontsize',8); %plot only 97
        plot(lon_st(i),lat_st(i),'.','color','k'); %plot only 97
    end
end
caxis([30 34])
colormap jet
saveas(fig,'2019_obs_bot_salt.png');

% nearest neigbor
% temp_intf_surf = temp_int_surf;
% salt_intf_surf = salt_int_surf;
% temp_intf_bot = temp_int_bot;
% salt_intf_bot = salt_int_bot;
% % temp_intf_surf(isnan(temp_int_surf))=interpextrap2(lon(~isnan(temp_int_surf)), lat(~isnan(temp_int_surf)),temp_int_surf(~isnan(temp_int_surf)),lon(isnan(temp_int_surf)),lon(isnan(temp_int_surf)));
% % temp_intf_surf(isnan(temp_int_surf))=griddata(lon(~isnan(temp_int_surf)), lat(~isnan(temp_int_surf)),temp_int_surf(~isnan(temp_int_surf)),lon(isnan(temp_int_surf)),lon(isnan(temp_int_surf)),'natural');
% salt_intf_surf=griddata(lon, lat, salt_int_surf,lon,lat,'nearest');
% temp_intf_bot=griddata(lon, lat,temp_int_bot,lon,lat,'nearest');
% salt_intf_bot=griddata(lon, lat,salt_int_bot,lon,lat,'nearest');

cd D:\RIST\RIST_±¤¾ç¸¸\vaildation\pic\t2_tem\horizontal

fig=figure; hold on;
pcolor(lon,lat,temp_int_surf.*(mask./mask)); shading flat; c = colorbar('Ticks',[0:3:15]);
% c.Label.String = '^oC';
c.Title.String = '^oC';
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
caxis([0 15]);
title(['surface temperature OBS'],'fontsize',30);set(gca,'FontSize',16);
xlabel('Lon','fontsize', 16);ylabel('Lat','fontsize',16);
for i=1:length(lat_obs)
%         %text(lon_obs(i)-0.001,lat_obs(i)-0.005, num2str(i),'color','k','fontsize',8); %plot only 97
        plot(lon_obs(i),lat_obs(i),'+','color','k','linew',2); %plot only 97
end
saveas(fig,'surf_obs_temp.png');

fig=figure; hold on;
pcolor(lon,lat,temp_int_bot.*(mask./mask)); shading flat; c = colorbar('Ticks',[0:3:15]);
% c.Label.String = '^oC';
c.Title.String = '^oC';
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
caxis([0 15]);
title(['bottom temperature OBS'],'fontsize',30);set(gca,'FontSize',16);
xlabel('Lon','fontsize', 16);ylabel('Lat','fontsize',16);
for i=1:length(lat_obs)
        %text(lon_obs(i)-0.001,lat_obs(i)-0.005, num2str(i),'color','k','fontsize',8); %plot only 97
        plot(lon_obs(i),lat_obs(i),'+','color','k','linew',2); %plot only 97
end
saveas(fig,'bot_obs_temp.png');

fig=figure; hold on;
pcolor(lon,lat,salt_int_surf.*(mask./mask)); shading flat;  c = colorbar('Ticks',[20:5:35]);
c.Title.String = 'ppt';
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
caxis([20 35]);
title(['surface salt OBS'],'fontsize',30);set(gca,'FontSize',16);
xlabel('Lon','fontsize', 16);ylabel('Lat','fontsize',16);
for i=1:length(lat_obs)
        %text(lon_obs(i)-0.001,lat_obs(i)-0.005, num2str(i),'color','k','fontsize',8); %plot only 97
        plot(lon_obs(i),lat_obs(i),'+','color','k','linew',2); %plot only 97
end
saveas(fig,'surf_obs_salt.png');

fig=figure; hold on;
pcolor(lon,lat,salt_int_bot.*(mask./mask)); shading flat; c = colorbar('Ticks',[20:5:35]);
c.Title.String = 'ppt';
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(lon)) 127.9])
caxis([20 35]);
title(['bottom salt OBS'],'fontsize',30);set(gca,'FontSize',16);
xlabel('Lon','fontsize', 16);ylabel('Lat','fontsize',16);
for i=1:length(lat_obs)
        %text(lon_obs(i)-0.001,lat_obs(i)-0.005, num2str(i),'color','k','fontsize',8); %plot only 97
        plot(lon_obs(i),lat_obs(i),'+','color','k','linew',2); %plot only 97
end
saveas(fig,'bot_obs_salt.png');




