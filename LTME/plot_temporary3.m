close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
obs = [ 2, 4, 6, 7, 8, 10, 11, 13];
for i = 1: length(obs)
name_tag_1(i) = obs(i); 
end

% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('1983_report_obs.xlsx','gy_1982','');

%%% column 1 : st., 2: depth, 3: time, 4: temp, 5: salt, 6: DO, 
%%% 7: ph, 8: NH4, 9: NO3

txt_matc_p = raw_p(1:end,1); % name list
txt_date_p = raw_p(1:end,3); % date list
txt_depth_p = raw_p(1:end,2); % depth list

% txt_cut=txt_p(3:end,:);

no_depth_nan_p = find(isnan(txt_depth_p)==0); % detect no NaN
depth_nan_p = find(isnan(txt_depth_p)==1); % detect NaN
for i = 1:length(name_tag_1)
   pick_st{i,1} = find(txt_matc_p(no_depth_nan_p) == name_tag_1(i));
   pick_bio_st{i,1} = find(txt_matc_p(depth_nan_p) == name_tag_1(i));
end

clearvars temp_p salt_p DO_p ph_p NH4_p NO3_p
temp_p = raw_p(:,4); 
salt_p = raw_p(:,5);
DO_p = raw_p(:,6);
ph_p = raw_p(:,7);
NH4_p = raw_p(:,8);
NO3_p = raw_p(:,9);
for i = 1:length(name_tag_1)
  time_obs{i,1} = txt_date_p(no_depth_nan_p(pick_st{i,1})); 
  depth_obs{i,1} = txt_depth_p(no_depth_nan_p(pick_st{i,1}));
  temp_obs{i,1} = temp_p(no_depth_nan_p(pick_st{i,1}));
  salt_obs{i,1} = salt_p(no_depth_nan_p(pick_st{i,1}));
  DO_obs{i,1} = DO_p(no_depth_nan_p(pick_st{i,1}));
  ph_obs{i,1} = ph_p(no_depth_nan_p(pick_st{i,1}));
  NH4_obs{i,1} = NH4_p(depth_nan_p(pick_bio_st{i,1}));
  NO3_obs{i,1} = NO3_p(depth_nan_p(pick_bio_st{i,1}));
end

clearvars diff_ind
for i = 1:length(name_tag_1)
    for j = 1:length(time_obs{i,1}(1:end))
         if j == 1
             clearvars temp_time
             temp_time = time_obs{i,1}(j);
             ij=0;
         end
         if temp_time ~= time_obs{i,1}(j);
            ij = ij+1;
            diff_ind{i}(ij) = j;
            clearvars temp_time
            temp_time = time_obs{i,1}(j);
         end
end
end


for i = 1:length(name_tag_1)
  for zx = 1:length(diff_ind{i})+1
      if zx == 1
          len_t = length(1:diff_ind{i}(zx)-1);
          t_obs{i,zx} = time_obs{i,1}(1:len_t);
          depth_obs_f{i,zx} = depth_obs{i,1}(1:len_t);
          temp_obs_f{i,zx} = temp_obs{i,1}(1:len_t);
          salt_obs_f{i,zx} = salt_obs{i,1}(1:len_t);
      else
          t_obs{i,zx} = time_obs{i,1}(len_t*(zx-1)+1:len_t*zx);
          depth_obs_f{i,zx} = depth_obs{i,1}(1:len_t);
          temp_obs_f{i,zx} = temp_obs{i,1}(1:len_t);
          salt_obs_f{i,zx} = salt_obs{i,1}(1:len_t);
 
      end
  end
%   DO_obs{i,1};
%   ph_obs{i,1};
%   NH4_obs{i,1};
%   NO3_obs{i,1};
end

% % obs_st
[raw_st txt_st]=xlsread('1983_report_obs.xlsx','gy_1982_st','');
lon_obs = raw_st(:,2); lat_obs = raw_st(:,3); 

% 1982, 5mth, 7mth, 10mth

% make 1980~present
k=0
for i = 1982:1982
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

clearvars em em_seanson k kk
k=[5, 7, 10, 12];
for i = 1:length(k)
    em{i} = sum(eom_d(1:k(i)-1))+1:sum(eom_d(1:k(i)));
end    

kk=[1:12];
for i = 1:length(kk)
    if i == 1
        em_season{i} = 1:sum(eom_d(1:kk(i)));
    else
        em_season{i} = sum(eom_d(1:kk(i)-1))+1:sum(eom_d(1:kk(i)));
    end
end

% model grid
lon = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_0001.nc'],'lon_rho');
lat = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_0001.nc'],'lat_rho');
h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_0001.nc'],'h');
mask =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_0001.nc'],'mask_rho');
% h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_0001.nc'],'h');

%% find nearest point on 

% % diff_lon=lon - lon_obs(1);
% diff_lon=[lon(1:end,1)] - lon_obs(1);
% min(abs(diff_lon))
% find(min(abs(diff_lon))==diff_lon)
% 
% % diff_lat=lat - lat_obs(1);
% diff_lat=[lat(1,1:end)] - lat_obs(1);
% min(abs(diff_lat))
% find(min(abs(diff_lat))==diff_lat)

for i = 1:length(lon_obs)
clearvars temp_2d
    temp_2d = sqrt((lon - lon_obs(i)).^2 + (lat - lat_obs(i)).^2);
    near_point(i) =find(min(min(temp_2d))==temp_2d);
    dist_2d(:,:,i) = temp_2d;
end

for i = 1:length(near_point)
    clearvars row col
    [row,col]=ind2sub(size(lon),near_point(i)); %% 1d to 2d index
    near_point_2d(i,:) = [row, col]; 
end

% min(min(dist_2d))
% find(min(min(dist_2d))==dist_2d)

for i = 1:length(k)
     for j = 1:length(em{i})
        temp(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em{i}(j),'%04d'),'.nc'],'temp');
        salt(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em{i}(j),'%04d'),'.nc'],'salt');
        zeta(:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em{i}(j),'%04d'),'.nc'],'zeta');
     end
end

for i = 1:12
     for j = 1:length(em_season{i}) 
        temp_all(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'temp');
        salt_all(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'salt');
        zeta_all(:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'zeta');
     end
end

% for i = 1:3
%      for j = 1:length(em_season{i}) 
%         no3(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NO3');
%         NH4(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NH4');
%         DO(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1982_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'oxygen');
%      end
% end


for i = 1:length(k)
        mth_temp(:,:,:,i) = sum(squeeze(temp(:,:,:,i,1:end)),4)./length(em{i});
        mth_salt(:,:,:,i) = sum(squeeze(salt(:,:,:,i,1:end)),4)./length(em{i});
        mth_zeta(:,:,i) = sum(squeeze(zeta(:,:,i,1:end)),3)./length(em{i});
end


for i = 1:12
        mth12_temp(:,:,:,i) = sum(squeeze(temp_all(:,:,:,i,1:end)),4)./length(em_season{i});
        mth12_salt(:,:,:,i) = sum(squeeze(salt_all(:,:,:,i,1:end)),4)./length(em_season{i});
        mth12_zeta(:,:,i) = sum(squeeze(zeta_all(:,:,i,1:end)),3)./length(em_season{i});
end


% for i = 1:12
%         mth_no3(:,:,:,i) = sum(squeeze(no3(:,:,:,i,1:end)),4)./length(em_season{i});
%         mth_NH4(:,:,:,i) = sum(squeeze(NH4(:,:,:,i,1:end)),4)./length(em_season{i});
%         mth_DO(:,:,:,i) = sum(squeeze(DO(:,:,:,i,1:end)),4)./length(em_season{i});
% end


figure; pcolor(lon,lat,squeeze(mth_temp(:,:,20,1))); shading flat;

clearvars z
theta_s = 1; theta_b=1;hc=1;N=20;
for i = 1:length(k)
    clearvars temp_z
    temp_z=zlevs(h.*(mask./mask), squeeze(mth_zeta(:,:,i)), theta_s, theta_b, hc, N, 'r');
    z(:,:,:,i)=temp_z;
end
% pcolor(lon,lat,squeeze(z(1,:,:))); shading flat; colorbar
% figure; pcolor(lon,lat,-h.*(mask./mask)-squeeze(z(1,:,:))); shading flat; colorbar

clearvars z_mod
for i =1:length(near_point)
    for j = 1:length(k)
        z_mod(:,i,j)=squeeze(z(:,near_point_2d(i,1),near_point_2d(i,2),j));
        temp_mod(:,i,j)=squeeze(mth_temp(near_point_2d(i,1),near_point_2d(i,2),:,j));
        salt_mod(:,i,j)=squeeze(mth_salt(near_point_2d(i,1),near_point_2d(i,2),:,j));
    end
end


clearvars z_12
theta_s = 1; theta_b=1;hc=1;N=20;
for i = 1:12
    clearvars temp_z
    temp_z=zlevs(h.*(mask./mask), squeeze(mth12_zeta(:,:,i)), theta_s, theta_b, hc, N, 'r');
    z_12(:,:,:,i)=temp_z;
end

clearvars z_mod
for i =1:length(near_point)
    for j = 1:12
        z_mod_12(:,i,j)=squeeze(z_12(:,near_point_2d(i,1),near_point_2d(i,2),j));
        temp_mod_12(:,i,j)=squeeze(mth12_temp(near_point_2d(i,1),near_point_2d(i,2),:,j));
        salt_mod_12(:,i,j)=squeeze(mth12_salt(near_point_2d(i,1),near_point_2d(i,2),:,j));
    end
end


return



for ff = 1:length(obs) %st.
    f=obs(ff);
    clearvars title_temp
    for ef = 1:4 % month
        clearvars t_mon 
        t_mon = k(ef);
        figure; plot(temp_obs_f{ff,ef}, -depth_obs_f{ff,ef},'color','r','linew',2); hold on;
        plot(temp_mod_12(:,f,t_mon),z_mod_12(:,f,t_mon),'b','linew',2);

    for i = 1:size(temp,5) % daily
        plot(squeeze(temp(near_point_2d(f,1),near_point_2d(f,2),:,ef,i)),z_mod_12(:,f,t_mon),'linew',2,'color',[169/255 169/255 169/255]);
    end
        plot(temp_mod_12(:,f,t_mon),z_mod_12(:,f,t_mon),'b','linew',2);
        plot(temp_obs_f{ff,ef}, -depth_obs_f{ff,ef},'*','color','r','linew',2);
        legend('location','southeast','OBS','model_m_o_n','model_d_a_y')
        grid on;
        xlim([5 30]);
        ylim([-32 0]);
        xlabel('temp(^oC)','fontsize',15)
        ylabel('depth(m)','fontsize',15)
        title_temp = ['temp compare st.', num2str(f),' - ' num2str(t_mon), 'mth'];
        title(title_temp)
        set(gca,'fontsize',15,'fontweight','bold')
        print('-dpng',[title_temp,'.png']);        
    end
   hold off;
end
close all


for ff = 1:length(obs) %st.
    f=obs(ff);
    clearvars title_temp
    for ef = 1:4 % month
        clearvars t_mon 
        t_mon = k(ef);
        figure; plot(salt_obs_f{ff,ef}, -depth_obs_f{ff,ef},'color','r','linew',2); hold on;
        plot(salt_mod_12(:,f,t_mon),z_mod_12(:,f,t_mon),'b','linew',2);

    for i = 1:size(salt,5) % daily
        plot(squeeze(salt(near_point_2d(f,1),near_point_2d(f,2),:,ef,i)),z_mod_12(:,f,t_mon),'linew',2,'color',[169/255 169/255 169/255]);
    end
        plot(salt_mod_12(:,f,t_mon),z_mod_12(:,f,t_mon),'b','linew',2);
        plot(salt_obs_f{ff,ef}, -depth_obs_f{ff,ef},'*','color','r','linew',2);
        legend('location','southwest','OBS','model_m_o_n','model_d_a_y')
        grid on;
        xlim([13 33.5]);
        ylim([-32 0]);
        xlabel('salt(g/kg)','fontsize',15)
        ylabel('depth(m)','fontsize',15)
        title_temp = ['salt compare st.', num2str(f),' - ' num2str(t_mon), 'mth'];
        title(title_temp)
        set(gca,'fontsize',15,'fontweight','bold')
        print('-dpng',[title_temp,'.png']);
    end
   hold off;
end
close all



% figure;
% for ff = 1:length(obs) %st.
%     f=obs(ff);
%     for ef = 1:4 % month
%         clearvars t_mon 
%         t_mon = k(ef);
%          plot(temp_obs_f{ff,ef}, -depth_obs_f{ff,ef},'r','linew',2); hold on;
%         plot(temp_mod_12(:,f,t_mon),z_mod_12(:,f,t_mon),'b','linew',2);
% 
%     for i = 1:size(temp,5) % daily
%         plot(squeeze(temp(near_point_2d(f,1),near_point_2d(f,2),:,ef,i)),z_mod_12(:,f,t_mon),'linew',2,'color',[169/255 169/255 169/255]);
%     end
%         plot(temp_mod_12(:,f,t_mon),z_mod_12(:,f,t_mon),'b','linew',2);
%         legend('OBS','model_m_o_n','model_d_a_y')
%         grid on;
%     end
% end
% 
% set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
%     'PaperUnits', 'Inches', 'PaperSize', [9, 7])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%% plot- DO vs. temps diff(surf-bot)


close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
obs = [ 18, 2];
for i = 1: length(obs)
name_tag_1(i) = obs(i); 
end

% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('1983_report_obs.xlsx','gy_1986','');

return
%%% column 1 : st., 2: depth, 3: time, 4: temp, 5: salt, 6: DO, 
%%% 7: ph, 8: NH4, 9: NO3

txt_matc_p = raw_p(1:end,1); % name list
txt_date_p = txt_p(2:end,1); % date list
% txt_depth_p = raw_p(1:end,2); % depth list

% txt_cut=txt_p(3:end,:);

no_tag_nan_p = find(isnan(txt_matc_p)==0); % detect no NaN
tag_nan_p = find(isnan(txt_matc_p)==1); % detect NaN
for i = 1:length(name_tag_1)
   pick_st{i,1} = find(txt_matc_p(no_tag_nan_p) == name_tag_1(i));
end
pick_st{2,1} = tag_nan_p;

for i = 1:length(name_tag_1)
obs_loc(i,1)=raw_p(pick_st{i,1}(1),2);
obs_loc(i,2)=raw_p(pick_st{i,1}(1),3);
end

clearvars temp_p salt_p DO_p ph_p NH4_p NO3_p
tt_date=char(txt_date_p(:,1));
temp_p = raw_p(:,4); 
salt_p = raw_p(:,5);
DO_p = raw_p(:,6);
ph_p = raw_p(:,8);
NH4_p = raw_p(:,11);
NO3_p = raw_p(:,10);
for i = 1:length(name_tag_1) 
  temp_obs{i,1} = temp_p((pick_st{i,1}));
  salt_obs{i,1} = salt_p((pick_st{i,1}));
  DO_obs{i,1} = DO_p((pick_st{i,1}));
  ph_obs{i,1} = ph_p((pick_st{i,1}));
  NH4_obs{i,1} = NH4_p((pick_st{i,1}));
  NO3_obs{i,1} = NO3_p((pick_st{i,1}));
end

clearvars time_obs
for i = 1:51
  time_obs{i} = {tt_date(pick_st{1,1}(i),6:end-3)};
end

%% time match (hr)
leapyear(1986)
% make 1980~present
k=0
for i = 1986:1986
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

clearvars em em_seanson k kk
k=[3];
for i = 1:length(k)
    em{i} = sum(eom_d(1:k(i)-1))+1:sum(eom_d(1:k(i)));
end    

% make hr
% 1986.03.11 09:00 ~ 03.12 10:00
% 1986.03.18 10:00 ~ 03.19 10:00

clearvars hr_em
% hr_em(1,:)=em{1}(10)*24 + 9 +(0:25) +9; % UTC to local time +9
% hr_em(1,end+1:end+length(0:24))=em{1}(17)*24 + 9 +(0:24) +9;
hr_em(1,:)=em{1}(10)*24 +(0:25) +9; % UTC to local time +9
hr_em(1,end+1:end+length(0:24))=em{1}(17)*24 +(0:24) +9;

% % obs_st
lon_obs = obs_loc(:,1); lat_obs = obs_loc(:,2); 

% 1982, 5mth, 7mth, 10mth



% model grid
lon = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'lon_rho');
lat = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'lat_rho');
h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'h');
mask =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'mask_rho');
% h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'h');

%% find nearest point on 

% % diff_lon=lon - lon_obs(1);
% diff_lon=[lon(1:end,1)] - lon_obs(1);
% min(abs(diff_lon))
% find(min(abs(diff_lon))==diff_lon)
% 
% % diff_lat=lat - lat_obs(1);
% diff_lat=[lat(1,1:end)] - lat_obs(1);
% min(abs(diff_lat))
% find(min(abs(diff_lat))==diff_lat)

for i = 1:length(lon_obs)
clearvars temp_2d
    temp_2d = sqrt((lon - lon_obs(i)).^2 + (lat - lat_obs(i)).^2);
    near_point(i) =find(min(min(temp_2d))==temp_2d);
    dist_2d(:,:,i) = temp_2d;
end

for i = 1:length(near_point)
    clearvars row col
    [row,col]=ind2sub(size(lon),near_point(i)); %% 1d to 2d index
    near_point_2d(i,:) = [row, col]; 
end

% min(min(dist_2d))
% find(min(min(dist_2d))==dist_2d)

for i = 1:length(hr_em)
        temp(:,:,:,i) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_his_',num2str(hr_em(i),'%04d'),'.nc'],'temp');
        salt(:,:,:,i) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_his_',num2str(hr_em(i),'%04d'),'.nc'],'salt');
        zeta(:,:,i) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_his_',num2str(hr_em(i),'%04d'),'.nc'],'zeta');
end


% for i = 1:3
%      for j = 1:length(em_season{i}) 
%         no3(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NO3');
%         NH4(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NH4');
%         DO(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'oxygen');
%      end
close all; clear; clc;   % -v3

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
obs = [ 18, 2];
for i = 1: length(obs)
name_tag_1(i) = obs(i); 
end

% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('1983_report_obs.xlsx','gy_1986','');

return
%%% column 1 : st., 2: depth, 3: time, 4: temp, 5: salt, 6: DO, 
%%% 7: ph, 8: NH4, 9: NO3

txt_matc_p = raw_p(1:end,1); % name list
txt_date_p = txt_p(2:end,1); % date list
% txt_depth_p = raw_p(1:end,2); % depth list

% txt_cut=txt_p(3:end,:);

no_tag_nan_p = find(isnan(txt_matc_p)==0); % detect no NaN
tag_nan_p = find(isnan(txt_matc_p)==1); % detect NaN
for i = 1:length(name_tag_1)
   pick_st{i,1} = find(txt_matc_p(no_tag_nan_p) == name_tag_1(i));
end
pick_st{2,1} = tag_nan_p;

for i = 1:length(name_tag_1)
obs_loc(i,1)=raw_p(pick_st{i,1}(1),2);
obs_loc(i,2)=raw_p(pick_st{i,1}(1),3);
end

clearvars temp_p salt_p DO_p ph_p NH4_p NO3_p
tt_date=char(txt_date_p(:,1));
temp_p = raw_p(:,4); 
salt_p = raw_p(:,5);
DO_p = raw_p(:,6);
ph_p = raw_p(:,8);
NH4_p = raw_p(:,11);
NO3_p = raw_p(:,10);
for i = 1:length(name_tag_1) 
  temp_obs{i,1} = temp_p((pick_st{i,1}));
  salt_obs{i,1} = salt_p((pick_st{i,1}));
  DO_obs{i,1} = DO_p((pick_st{i,1}));
  ph_obs{i,1} = ph_p((pick_st{i,1}));
  NH4_obs{i,1} = NH4_p((pick_st{i,1}));
  NO3_obs{i,1} = NO3_p((pick_st{i,1}));
end

clearvars time_obs
for i = 1:51
  time_obs{i} = {tt_date(pick_st{1,1}(i),6:end-3)};
end

%% time match (hr)
leapyear(1986)
% make 1980~present
k=0
for i = 1986:1986
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

clearvars em em_seanson k kk
k=[3];
for i = 1:length(k)
    em{i} = sum(eom_d(1:k(i)-1))+1:sum(eom_d(1:k(i)));
end    

% make hr
% 1986.03.11 09:00 ~ 03.12 10:00
% 1986.03.18 10:00 ~ 03.19 10:00

clearvars hr_em
% hr_em(1,:)=em{1}(10)*24 + 9 +(0:25) +9; % UTC to local time +9
% hr_em(1,end+1:end+length(0:24))=em{1}(17)*24 + 9 +(0:24) +9;
hr_em(1,:)=em{1}(10)*24 +(0:25) +9; % UTC to local time +9
hr_em(1,end+1:end+length(0:24))=em{1}(17)*24 +(0:24) +9;

% % obs_st
lon_obs = obs_loc(:,1); lat_obs = obs_loc(:,2); 

% 1982, 5mth, 7mth, 10mth



% model grid
lon = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'lon_rho');
lat = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'lat_rho');
h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'h');
mask =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'mask_rho');
% h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'h');

%% find nearest point on 

% % diff_lon=lon - lon_obs(1);
% diff_lon=[lon(1:end,1)] - lon_obs(1);
% min(abs(diff_lon))
% find(min(abs(diff_lon))==diff_lon)
% 
% % diff_lat=lat - lat_obs(1);
% diff_lat=[lat(1,1:end)] - lat_obs(1);
% min(abs(diff_lat))
% find(min(abs(diff_lat))==diff_lat)

for i = 1:length(lon_obs)
clearvars temp_2d
    temp_2d = sqrt((lon - lon_obs(i)).^2 + (lat - lat_obs(i)).^2);
    near_point(i) =find(min(min(temp_2d))==temp_2d);
    dist_2d(:,:,i) = temp_2d;
end

for i = 1:length(near_point)
    clearvars row col
    [row,col]=ind2sub(size(lon),near_point(i)); %% 1d to 2d index
    near_point_2d(i,:) = [row, col]; 
end

% min(min(dist_2d))
% find(min(min(dist_2d))==dist_2d)

for i = 1:length(hr_em)
        temp(:,:,:,i) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_his_',num2str(hr_em(i),'%04d'),'.nc'],'temp');
        salt(:,:,:,i) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_his_',num2str(hr_em(i),'%04d'),'.nc'],'salt');
        zeta(:,:,i) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_his_',num2str(hr_em(i),'%04d'),'.nc'],'zeta');
end


% for i = 1:3
%      for j = 1:length(em_season{i}) 
%         no3(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NO3');
%         NH4(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NH4');
%         DO(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'oxygen');
%      end
% end

for i =1:length(near_point)
        temp_mod(i,:)=squeeze(temp(near_point_2d(i,1),near_point_2d(i,2),20,:));
        salt_mod(i,:)=squeeze(salt(near_point_2d(i,1),near_point_2d(i,2),20,:));
end

% figure;
% for i = 1:1 %st.
%   plot(temp_obs{i},'color','r','linew',2); hold on;
%   plot(temp_mod(i,:)-0.7152,'b','linew',2);
% end

figure;
for i = 1:1 %st.
  plot(temp_obs{i},'color','r','linew',2); hold on;
  plot(temp_mod(i,:),'b','linew',2);
  legend('location','southeast','OBS_h_r_1_8','model_h_r_1_8')
  grid on;
        xlim([1 inf]);
        set(gca,'xtick',[1 27 51])
        set(gca,'xticklabels',{'3/11 09','3/18 10','3/19 10'})
        ylim([0 10]);
        xlabel('time(mmddhh)','fontsize',15)
        ylabel('temp(^oC)','fontsize',15)
        title_temp = ['temp compare st.', num2str(obs(i)),' - 1986'];
        title(title_temp)
        set(gca,'fontsize',15,'fontweight','bold')
        print('-dpng',[title_temp,'.png']);
        ax=gca;
        ax.XAxis.MinorTick='on'
        ax.XAxis.MinorTickValues=1:1:51
        set(ax,'xminorgrid','on')
end

line_spc = {'-','--'};
figure;
for i = 1:2 %st.
  plot(salt_obs{i},line_spc{i},'color','r','linew',2); hold on;
  plot(salt_mod(i,:),line_spc{i},'color','b','linew',2);
  legend('location','southwest','OBS_h_r_1_8','model_h_r_1_8','OBS_h_r_2_A','model_h_r_2_A')
  grid on;
        xlim([1 inf]);
        set(gca,'xtick',[1 27 51])
        set(gca,'xticklabels',{'3/11 09','3/18 10','3/19 10'})
        xlabel('time(mmddhh)','fontsize',15)
        ylabel('salt(g/kg)','fontsize',15)
        title_temp = ['salt compare st. 18 & 2A - 1986'];
        title(title_temp)
        set(gca,'fontsize',15,'fontweight','bold')
        grid on;
        ax=gca;
        ax.XAxis.MinorTick='on'
        ax.XAxis.MinorTickValues=1:1:51
        set(ax,'xminorgrid','on')
        print('-dpng',[title_temp,'.png']);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear; clc;   % -v3 1986~1987

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 

% port
obs = [1:5];
for i = 1: length(obs)
name_tag_1(i) = obs(i); 
end

% combining the tag and outter point excluding
% name_tag = name_tag_1'; 
% name_tag{end+1:end+length(name_tag_2)} = name_tag_2; 
% name_tag{end+1:end+length(name_tag_3)} = name_tag_3; 
% size_tag = length(name_tag);


%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('1983_report_obs.xlsx','gy_1986_1987','');

%%% column 1 : st., 2: depth, 3: time, 4: temp, 5: salt, 6: DO, 
%%% 7: ph, 8: NH4, 9: NO3

txt_matc_p = raw_p(1:end,2); % name list
txt_date_p = {'1986.09','1986.11','1987.02','1987.04','1987.06'}; % date list
% txt_depth_p = raw_p(1:end,2); % depth list

% txt_cut=txt_p(3:end,:);

no_tag_nan_p = find(isnan(txt_matc_p)==0); % detect no NaN
tag_nan_p = find(isnan(txt_matc_p)==1); % detect NaN
for i = 1:length(name_tag_1)
   pick_st{i,1} = find(txt_matc_p(no_tag_nan_p) == name_tag_1(i));
end


for i = 1:length(name_tag_1)
obs_loc(i,1)=raw_p(pick_st{i,1}(1),3);
obs_loc(i,2)=raw_p(pick_st{i,1}(1),4);
end

clearvars temp_p salt_p DO_p chla_p NO3_p
temp_p = raw_p(:,5); 
salt_p = raw_p(:,6);
DO_p = raw_p(:,7);
chla_p = raw_p(:,8);
NO3_p = raw_p(:,9);
for i = 1:length(name_tag_1) 
  temp_obs{i,1} = temp_p((pick_st{i,1}));
  salt_obs{i,1} = salt_p((pick_st{i,1}));
  DO_obs{i,1} = DO_p((pick_st{i,1}));
  NO3_obs{i,1} = NO3_p((pick_st{i,1}));
  chla_obs{i,1} = chla_p((pick_st{i,1}));
end


%% time match 
leapyear(1986:1987)
% make 1980~present
k=0
for i = 1986:1986
            k=k+1;
    for j = 1:12
        eom_d_f(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0
for i = 1987:1987
            k=k+1;
    for j = 1:12
        eom_d_b(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end


clearvars em em_seanson k_f kk
k_f=[9 11];
for i = 1:length(k_f)
    em_f{i} = sum(eom_d_f(1:k_f(i)-1))+1:sum(eom_d_f(1:k_f(i)));
end

clearvars em em_seanson k_b kk
k_b=[2 4 6];
for i = 1:length(k_b)
    em_b{i} = sum(eom_d_b(1:k_b(i)-1))+1:sum(eom_d_b(1:k_b(i)));
end    


% % obs_st
lon_obs = obs_loc(:,1); lat_obs = obs_loc(:,2); 

% 1982, 5mth, 7mth, 10mth



% model grid
lon = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'lon_rho');
lat = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'lat_rho');
h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'h');
mask =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'mask_rho');
% h =  ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_0001.nc'],'h');

%% find nearest point on 

% % diff_lon=lon - lon_obs(1);
% diff_lon=[lon(1:end,1)] - lon_obs(1);
% min(abs(diff_lon))
% find(min(abs(diff_lon))==diff_lon)
% 
% % diff_lat=lat - lat_obs(1);
% diff_lat=[lat(1,1:end)] - lat_obs(1);
% min(abs(diff_lat))
% find(min(abs(diff_lat))==diff_lat)
clearvars near_point*
for i = 1:length(lon_obs)
clearvars temp_2d 
    temp_2d = sqrt((lon - lon_obs(i)).^2 + (lat - lat_obs(i)).^2);
    near_point{i} =find(min(min(temp_2d))==temp_2d);
    dist_2d(:,:,i) = temp_2d;
end

for i = 1:length(near_point)
    clearvars row col
    [row,col]=ind2sub(size(lon),near_point{i}(1)); %% 1d to 2d index
    near_point_2d(i,:) = [row, col]; 
end

% min(min(dist_2d))
% find(min(min(dist_2d))==dist_2d)

for i = 1:length(k_f)
     for j = 1:length(em_f{i})
        temp_f(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_f{i}(j),'%04d'),'.nc'],'temp');
        salt_f(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_f{i}(j),'%04d'),'.nc'],'salt');
        zeta_f(:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_f{i}(j),'%04d'),'.nc'],'zeta');
     end
end

for i = 1:length(k_b)
     for j = 1:length(em_b{i})
        temp_b(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1987_regrid_v3/ocean_avg_',num2str(em_b{i}(j),'%04d'),'.nc'],'temp');
        salt_b(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1987_regrid_v3/ocean_avg_',num2str(em_b{i}(j),'%04d'),'.nc'],'salt');
        zeta_b(:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1987_regrid_v3/ocean_avg_',num2str(em_b{i}(j),'%04d'),'.nc'],'zeta');
     end
end



% for i = 1:3
%      for j = 1:length(em_season{i}) 
%         no3(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NO3');
%         NH4(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'NH4');
%         DO(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/1986_regrid_v3/ocean_avg_',num2str(em_season{i}(j),'%04d'),'.nc'],'oxygen');
%      end

for i = 1:length(k_f)
        mth_temp_f(:,:,:,i) = sum(squeeze(temp_f(:,:,:,i,1:end)),4)./length(em_f{i});
        mth_salt_f(:,:,:,i) = sum(squeeze(salt_f(:,:,:,i,1:end)),4)./length(em_f{i});
        mth_zeta_f(:,:,i) = sum(squeeze(zeta_f(:,:,i,1:end)),3)./length(em_f{i});
end

for i = 1:length(k_b)
        mth_temp_b(:,:,:,i) = sum(squeeze(temp_b(:,:,:,i,1:end)),4)./length(em_b{i});
        mth_salt_b(:,:,:,i) = sum(squeeze(salt_b(:,:,:,i,1:end)),4)./length(em_b{i});
        mth_zeta_b(:,:,i) = sum(squeeze(zeta_b(:,:,i,1:end)),3)./length(em_b{i});
end

% merge the data (daily)
si_leng = size(temp_f);
day_temp = NaN(si_leng(1),si_leng(2),si_leng(3),length(k_f)+length(k_b),30);
day_salt = NaN(si_leng(1),si_leng(2),si_leng(3),length(k_f)+length(k_b),30);
day_zeta = NaN(si_leng(1),si_leng(2),length(k_f)+length(k_b),30);

day_temp(:,:,:,1:length(k_f),:) = temp_f; day_temp(:,:,:,length(k_f)+1:length(k_f)+length(k_b),:) = temp_b;
day_salt(:,:,:,1:length(k_f),:) = salt_f; day_salt(:,:,:,length(k_f)+1:length(k_f)+length(k_b),:) = salt_b;
day_zeta(:,:,1:length(k_f)) = mth_zeta_f; day_zeta(:,:,length(k_f)+1:length(k_f)+length(k_b)) = zeta_b;

% merge the data (monthly)
si_leng = size(mth_temp_f);
mth_temp = NaN(si_leng(1),si_leng(2),si_leng(3),length(k_f)+length(k_b));
mth_salt = NaN(si_leng(1),si_leng(2),si_leng(3),length(k_f)+length(k_b));
mth_zeta = NaN(si_leng(1),si_leng(2),length(k_f)+length(k_b));

mth_temp(:,:,:,1:length(k_f)) = mth_temp_f; mth_temp(:,:,:,length(k_f)+1:length(k_f)+length(k_b)) = mth_temp_b;
mth_salt(:,:,:,1:length(k_f)) = mth_salt_f; mth_salt(:,:,:,length(k_f)+1:length(k_f)+length(k_b)) = mth_salt_b;
mth_zeta(:,:,1:length(k_f)) = mth_zeta_f; mth_zeta(:,:,length(k_f)+1:length(k_f)+length(k_b)) = mth_zeta_b;

clearvars z_mer
theta_s = 1; theta_b=1;hc=1;N=20;
for i = 1:length(k_f)+length(k_b)
    clearvars temp_z
    temp_z=zlevs(h.*(mask./mask), squeeze(mth_zeta(:,:,i)), theta_s, theta_b, hc, N, 'r');
    z_mer(:,:,:,i)=temp_z;
end

clearvars temp_mod* salt_mod* z_mod*
% *_mod (N,st., time);
for i =1:length(near_point)
    for j = 1:length(k_f)+length(k_b) %time ind.
        temp_mod(:,i,j)=squeeze(mth_temp(near_point_2d(i,1),near_point_2d(i,2),:,j));
        salt_mod(:,i,j)=squeeze(mth_salt(near_point_2d(i,1),near_point_2d(i,2),:,j));
        z_mod(:,i,j)=squeeze(z_mer(:,near_point_2d(i,1),near_point_2d(i,2),j));
        temp_mod_raw(:,i,j,:)=squeeze(day_temp(near_point_2d(i,1),near_point_2d(i,2),:,j,:));
        salt_mod_raw(:,i,j,:)=squeeze(day_salt(near_point_2d(i,1),near_point_2d(i,2),:,j,:));
    end
end


% figure;
% for i = 1:1 %st.
%   plot(temp_obs{i},'color','r','linew',2); hold on;
%   plot(temp_mod(i,:)-0.7152,'b','linew',2);
% end
x_in = ones(5,30);
x_in(1,:) = 1; x_in(2,:) = 2; x_in(3,:) = 3; x_in(4,:) = 4; x_in(5,:) = 5;



for i = 1:length(near_point) %st.
figure;
%     for j = 1:length(k_f)+length(k_b) %time.
  plot(temp_obs{i},'color','r','linew',2); hold on;
  plot((squeeze(temp_mod(end,i,:)) + squeeze(temp_mod(11,i,:)) + squeeze(temp_mod(3,i,:)))./3,'b','linew',2); %surf
%   plot(squeeze(temp_mod(end,i,:)),'linew',2,'color',[169/255 169/255 169/255]); %surf
%   plot(squeeze(temp_mod(11,i,:)),'linew',2,'color',[169/255 169/255 169/255]); %mid
%   plot(squeeze(temp_mod(3,i,:)),'linew',2,'color',[169/255 169/255 169/255]); %bot
  plot(x_in, (squeeze(temp_mod_raw(end,i,:,:)) + squeeze(temp_mod_raw(11,i,:,:)) + squeeze(temp_mod_raw(3,i,:,:)))./3,'.','linew',2,'color','k'); %surf
  legend('location','north','OBS','model_m_o_n_t_h_l_y','model_d_a_i_l_y')
  grid on;
        xlim([1 inf]);
        set(gca,'xtick',[1 2 3 4 5])
        set(gca,'xticklabels',{'86.09','86.11','87.02','87.04','87.06'})
        ylim([0 28]);
        xlabel('time(yy.mm)','fontsize',15)
        ylabel('temp(^oC)','fontsize',15)
        title_temp = ['temp compare st.', num2str(obs(i)),' - 1986 ~ 1987'];
        title(title_temp)
        set(gca,'fontsize',15,'fontweight','bold')
        print('-dpng',[title_temp,'.png']);
        ax=gca;
        ax.XAxis.MinorTick='on'
        ax.XAxis.MinorTickValues=1:1:51
        set(ax,'xminorgrid','on')
        hold off;
end

for i = 1:length(near_point) %st.
figure;
%     for j = 1:length(k_f)+length(k_b) %time.
  plot(salt_obs{i},'color','r','linew',2); hold on;
  plot((squeeze(salt_mod(end,i,:)) + squeeze(salt_mod(11,i,:)) + squeeze(salt_mod(3,i,:)))./3,'b','linew',2); %surf
%   plot(squeeze(salt_mod(end,i,:)),'linew',2,'color',[169/255 169/255 169/255]); %surf
%   plot(squeeze(salt_mod(11,i,:)),'linew',2,'color',[169/255 169/255 169/255]); %mid
%   plot(squeeze(salt_mod(3,i,:)),'linew',2,'color',[169/255 169/255 169/255]); %bot
  plot(x_in, (squeeze(salt_mod_raw(end,i,:,:)) + squeeze(salt_mod_raw(11,i,:,:)) + squeeze(salt_mod_raw(3,i,:,:)))./3,'.','linew',2,'color','k'); %surf
  legend('location','south','OBS','model_m_o_n_t_h_l_y','model_d_a_i_l_y')
  grid on;
        xlim([1 inf]);
        set(gca,'xtick',[1 2 3 4 5])
        set(gca,'xticklabels',{'86.09','86.11','87.02','87.04','87.06'})
        ylim([20 35]);
        xlabel('time(yy.mm)','fontsize',15)
        ylabel('salt(PSU)','fontsize',15)
        title_salt = ['salt compare st.', num2str(obs(i)),' - 1986 ~ 1987'];
        title(title_salt)
        set(gca,'fontsize',15,'fontweight','bold')
        print('-dpng',[title_salt,'.png']);
        ax=gca;
        ax.XAxis.MinorTick='on'
        ax.XAxis.MinorTickValues=1:1:51
        set(ax,'xminorgrid','on')
        hold off;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
%% Gure river compare
% /data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY

% model grid
temp = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/river_1982_realts_biofennel_Gwangyang.nc'],'river_temp');
DO = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/river_1982_realts_biofennel_Gwangyang.nc'],'river_Oxyg');
NO3 = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/river_1982_realts_biofennel_Gwangyang.nc'],'river_NO3');
NH4 = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/input_GY/river_1982_realts_biofennel_Gwangyang.nc'],'river_NH4');

r=1;  % songjung river flow index

%% time match 
leapyear(1982)
% make 1980~present
k=0
for i = 1982:1982
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end


clearvars em em_seanson k_f kk
k_f=[9 10];
for i = 1:length(k_f)
    em{i} = sum(eom_d(1:k_f(i)-1))+1:sum(eom_d(1:k_f(i)));
end

for i =1:2
    if i == 1
        em_fix{i} = em{i}(12:end);
    else
        em_fix{i} = em{i}(1:12); 
    end
end

for i = 1:2
    temp_mod{i} = squeeze(temp(r,20,em_fix{i}));
end

%merge data
temp_mod_mer(1:length(temp_mod{1})) = temp_mod{1};
temp_mod_mer(length(temp_mod{1})+1:length(temp_mod{1})+length(temp_mod{2})) = temp_mod{2};

% clearvars his_ind
% his_ind(1:2,1) = [(em_fix{1}(1) - 1) * 24; (em_fix{2}(end) - 1) * 24 + 23;];

temp_mod_83 = (mean(temp_mod{1}) + mean(temp_mod{2})) ./2 ;


%% pick the row on the excel which has same name with tag
[raw_p txt_p]=xlsread('1983_report_obs.xlsx','sumjin','');

temp_obs = raw_p(1);
DO_obs = raw_p(3);
NH3_obs = raw_p(4);
NO3_obs = raw_p(5);


figure;
% h_c = plot([1:31],[ones(length(1:31)).*temp_obs],'color','r', ...
%     [1:31],[ones(length(1:31)).*temp_mod_83],'color','b','linew',2);
%     [1:31],temp_mod_mer,'color',[169/255 169/255 169/255],'linew',2);
plot([1:31],[ones(1,31).*temp_obs],'color','r','linew',2); hold on;
plot([1:31],[ones(1,31).*temp_mod_83],'color','b','linew',2);
%   plot(1:31,ones(length(1:31)).*temp_obs,'color','r','linew',2)
%   plot(1:31,ones(length(1:31)).*temp_mod_83,'color','b','linew',2)
  plot([1:31],temp_mod_mer,'linew',2,'color',[169/255 169/255 169/255]);
  legend('location','southwest','OBS','model_m_e_a_n','model_d_a_i_l_y')
  grid on;
        xlim([1 inf]);
        set(gca,'xtick',[1:5:31])
        set(gca,'xticklabels',{'09/12','09/17','09/22','09/27','10/02','10/07','10/12'})
        ylim([18 inf]);
        xlabel('time(yy/mm)','fontsize',15)
        ylabel('temp(^oC)','fontsize',15)
        title_temp = ['temp compare st.4 - 1982'];
        title(title_temp)
        set(gca,'fontsize',15,'fontweight','bold')
        print('-dpng',[title_temp,'.png']);
%         ax=gca;
%         ax.XAxis.MinorTick='on'
%         ax.XAxis.MinorTickValues=1:1:51
%         set(ax,'xminorgrid','on')
        hold off;

        

        

