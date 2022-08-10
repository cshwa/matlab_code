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
k=[5, 7, 10];
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

for i = 1:3
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


for i = 1:3
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
for i = 1:3
    clearvars temp_z
    temp_z=zlevs(h.*(mask./mask), squeeze(mth_zeta(:,:,i)), theta_s, theta_b, hc, N, 'r');
    z(:,:,:,i)=temp_z;
end
% pcolor(lon,lat,squeeze(z(1,:,:))); shading flat; colorbar
% figure; pcolor(lon,lat,-h.*(mask./mask)-squeeze(z(1,:,:))); shading flat; colorbar

clearvars z_mod
for i =1:length(near_point)
    for j = 1:3
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


figure; plot(temp_obs{1}(1:5), -depth_obs{1}(1:5),'r','linew',2); hold on;
plot(temp_mod_12(:,1,5),z_mod_12(:,1,5),'b','linew',2);
kj=1
for i = 1:size(temp,5)
    plot(squeeze(temp(near_point_2d(kj,1),near_point_2d(kj,2),:,1,i)),z_mod_12(:,1,5),'linew',2,'color',[169/255 169/255 169/255]);
end
plot(temp_mod_12(:,1,5),z_mod_12(:,1,5),'b','linew',2);
legend('OBS','model_m_o_n','model_d_a_y')
grid on;





% plot(mean(squeeze(temp_mod_12(:,1,3:5)),2), mean(squeeze(z_mod_12(:,1,3:5)),2),'k');


figure;
plot(temp_mod(:,1,2),z_mod(:,1,2)); hold on; plot(temp_obs{1}(6:10), -depth_obs{1}(6:10),'r')

figure; plot(temp_mod(:,1,3),z_mod(:,1,3)); hold on; plot(temp_obs{1}(11:15), -depth_obs{1}(11:15),'r')




%% make date to be 'yymm' form
for i = 1:length(name_tag_1)
clearvars temp temp_ymd temp_ymd_c temp_yymmdd temp_yymmdd_c
temp = char(date_st{i});
temp_ymd=temp(:,1:end-12);
temp_yymmdd=temp(:,1:end-9);

for j = 1:length(temp_ymd)
    temp_ymd_c{j} = temp_ymd(j,:);
    temp_yymmdd_c{j} = temp_yymmdd(j,:);
                    
end
date_ymd{i,1} = temp_ymd_c;
date_yymmdd{i,1} = temp_yymmdd_c;
end



k=0; m=0;
for i = 1:40
    l=0
        ref_yy(i,:)=[num2str(i+1979)];
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+1979) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+1979) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
    end
    end
end


% matched date 'yymm' form on surf & bottom
for j = 1:length(name_tag_1) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_ymd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
       end
    end 
end

indx_date_s = cell(1,length(ref_yymmdd));
indx_date_b = cell(1,length(ref_yymmdd));
for j = 1:length(name_tag_1) % st. axis
    for i = 1:length(ref_yymmdd) % date axis
       if  sum(strcmp(ref_yymmdd{i}, date_yymmdd{j})) ~= 0
           clearvars indx_date
           indx_date = find([strcmp(ref_yymmdd{i}, date_yymmdd{j})] == 1); 
           indx_date_s{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == 0));
           indx_date_b{j,i} = indx_date(find(depth_kodc{j}(1,indx_date) == max(depth_kodc{j}(1,indx_date))));
%        elseif sum(strcmp(ref_yymmdd{i}, date_yymmdd{j})) == 0
%            indx_date_s{j,i} = NaN;
%            indx_date_b{j,i} = NaN;
       end
    end
end


% clearvars indx_date_*
% for j = 1:3 % st. axis
%     for i = 1:length(ref_yymm) % date axis
%        if  sum(strcmp(ref_yymm{i}, date_ymd{j})) ~= 0
%            clearvars indx_date
%            indx_date = find([strcmp(ref_yymm{i}, date_ymd{j})] == 1);
%            if sum(indx_date(find(depth_kodc{j}(1,indx_date)==0))) ~= 0
%             indx_date_s(j,i) = indx_date(find(depth_kodc{j}(1,indx_date)==0));
%            elseif  sum(indx_date(find(depth_kodc{j}(1,indx_date)==0))) == 0
%             indx_date_s(j,i) = NaN;
%            end
%            
% %            indx_date_b(j,i) = indx_date(find(depth_kodc{j}(1,indx_date)==max(depth_kodc{j}(1,indx_date))));
%        end
%     end
% end


% make climate

%temp
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = temp_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        temp_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        temp_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = temp_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        temp_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        temp_bot(i,j) = NaN;
    end
    end
end

%salt
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = salt_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        salt_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        salt_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = salt_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        salt_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        salt_bot(i,j) = NaN;
    end
    end
end

%do
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = do_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        do_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        do_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = do_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        do_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        do_bot(i,j) = NaN;
    end
    end
end


%no3
for i = 1:size(indx_date_s,1)
    clearvars temp_kodc_k
    temp_kodc_k = no3_kodc{i};
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_kodc_k(indx_date_s{i,j}))) ~= 0     
        no3_sur(i,j) = nanmean(temp_kodc_k(indx_date_s{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_s{i,j}))) == 0 
        no3_sur(i,j) = NaN;
    end
    end
end

for i = 1:size(indx_date_b,1)
    clearvars temp_kodc_k
    temp_kodc_k = no3_kodc{i};
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_kodc_k(indx_date_b{i,j}))) ~= 0     
        no3_bot(i,j) = nanmean(temp_kodc_k(indx_date_b{i,j}));
    elseif sum(size(temp_kodc_k(indx_date_b{i,j}))) == 0 
        no3_bot(i,j) = NaN;
    end
    end
end

% check when they were overapped
%surface
for i = 1:size(indx_date_s,1)
    for j = 1:size(indx_date_s,2)
    if sum(size(temp_p(indx_date_s{i,j}))) > 2     
       disp(['i val = ', num2str(i)]); disp(['j val = ', num2str(j)]);
    end
    end
end

%bottom
for i = 1:size(indx_date_b,1)
    for j = 1:size(indx_date_b,2)
    if sum(size(temp_p(indx_date_b{i,j}))) > 2     
       disp(['i val = ', num2str(i)]); disp(['j val = ', num2str(j)]);
    end
    end
end

% save('KODC_data_monthly_20501.mat','*_bot','*_sur','ref_yymm', 'date_ymd');


%% interp >> to do regime shift test
% load KODC_data_monthly_20501.mat

do_sur = do_sur(1,:)';
do_bot = do_bot(1,:)';
no3_sur =  no3_sur(1,:)'.* 14;  %umol/L -> ug/L
no3_bot = no3_bot(1,:)' .* 14;
temp_sur = temp_sur(1,:)';
temp_bot = temp_bot(1,:)';
salt_sur = salt_sur(1,:)';
salt_bot = salt_bot(1,:)';

%% extract over 3sig
sig = 3; %% sigma
% sur
clearvars idx1
idx1 = find(isnan(do_sur) == 0);
regime_do=do_sur;
regime_do(find(do_sur > mean(do_sur(idx1)) + sig*std(do_sur(idx1))))=NaN;
regime_do(find(do_sur < mean(do_sur(idx1)) - sig*std(do_sur(idx1))))=NaN;
regm_do = nanmean(regime_do)

clearvars idx1
idx1 = find(isnan(no3_sur) == 0);
regime_no3=no3_sur;
regime_no3(find(no3_sur > mean(no3_sur(idx1)) + sig*std(no3_sur(idx1))))=NaN;
regime_no3(find(no3_sur < mean(no3_sur(idx1)) - sig*std(no3_sur(idx1))))=NaN;
regm_no3 = nanmean(regime_no3)


clearvars idx1
idx1 = find(isnan(temp_sur) == 0);
regime_temp=temp_sur;
regime_temp(find(temp_sur > mean(temp_sur(idx1)) + sig*std(temp_sur(idx1))))=NaN;
regime_temp(find(temp_sur < mean(temp_sur(idx1)) - sig*std(temp_sur(idx1))))=NaN;
regm_temp = nanmean(regime_temp)


clearvars idx1
idx1 = find(isnan(salt_sur) == 0);
regime_salt=salt_sur;
regime_salt(find(salt_sur > mean(salt_sur(idx1)) + sig*std(salt_sur(idx1))))=NaN;
regime_salt(find(salt_sur < mean(salt_sur(idx1)) - sig*std(salt_sur(idx1))))=NaN;
regm_salt = nanmean(regime_salt)

%bot
clearvars idx1
idx1 = find(isnan(do_bot) == 0);
regime_do_b=do_bot;
regime_do_b(find(do_bot > mean(do_bot(idx1)) + sig*std(do_bot(idx1))))=NaN;
regime_do_b(find(do_bot < mean(do_bot(idx1)) - sig*std(do_bot(idx1))))=NaN;
regm_do_b = nanmean(regime_do_b)

clearvars idx1
idx1 = find(isnan(no3_bot) == 0);
regime_no3_b=no3_bot;
regime_no3_b(find(no3_bot > mean(no3_bot(idx1)) + sig*std(no3_bot(idx1))))=NaN;
regime_no3_b(find(no3_bot < mean(no3_bot(idx1)) - sig*std(no3_bot(idx1))))=NaN;
regm_no3_b = nanmean(regime_no3_b)


clearvars idx1
idx1 = find(isnan(temp_bot) == 0);
regime_temp_b=temp_bot;
regime_temp_b(find(temp_bot > mean(temp_bot(idx1)) + sig*std(temp_bot(idx1))))=NaN;
regime_temp_b(find(temp_bot < mean(temp_bot(idx1)) - sig*std(temp_bot(idx1))))=NaN;
regm_temp_b = nanmean(regime_temp)


clearvars idx1
idx1 = find(isnan(salt_bot) == 0);
regime_salt_b=salt_bot;
regime_salt_b(find(salt_bot > mean(salt_bot(idx1)) + sig*std(salt_bot(idx1))))=NaN;
regime_salt_b(find(salt_bot < mean(salt_bot(idx1)) - sig*std(salt_bot(idx1))))=NaN;
regm_salt_b = nanmean(regime_salt)

%% extract regime mean
regime_do = regime_do - regm_do;
regime_no3 = regime_no3 - regm_no3;
regime_temp = regime_temp - regm_temp;
regime_salt = regime_salt - regm_salt;

regime_do_b = regime_do_b - regm_do_b;
regime_no3_b = regime_no3_b - regm_no3_b;
regime_temp_b = regime_temp_b - regm_temp_b;
regime_salt_b = regime_salt_b - regm_salt_b;

%% matched 366

% make 366 mm-dd
for i = 1:12
  eom_d(i) = eomday(1980,i); % 1980 is leap-yr
end

k=0
for i = 1:12
    l=0
    for j = 1:eom_d(i)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mmdd(k,:)=[num2str(i,'%02d') '-'  num2str(l,'%02d')]
    end
end

% make it to cell-array
for i = 1:length(mmdd)
    mmdd_c{i,1} = mmdd(i,:); % delete year
end

% pick matched date from water temp date
for i = 1:length(mmdd_c)
       mmdd_indx{i} = find([strcmp(mmdd_c{i}, ref_mmdd)] == 1)
end


%confirm how may data on each the days
size_c = [];
for i = 1:366; size_c(i)=length(mmdd_indx{i}); end
bar(size_c)

% make quasi_climate 1980~2019
for i = 1:length(mmdd_indx)
    if size(mmdd_indx{i},1) == 0     
        reg_clim_do(i) = NaN;
        reg_clim_no3(i) = NaN;
        reg_clim_temp(i) = NaN; 
        reg_clim_salt(i) = NaN;
        reg_clim_do_b(i) = NaN;
        reg_clim_no3_b(i) = NaN;
        reg_clim_temp_b(i) = NaN; 
        reg_clim_salt_b(i) = NaN;
    else
        reg_clim_do(i) = nanmean(regime_do(mmdd_indx{i}));
        reg_clim_no3(i) = nanmean(regime_no3(mmdd_indx{i}));
        reg_clim_temp(i) = nanmean(regime_temp(mmdd_indx{i}));
        reg_clim_salt(i) = nanmean(regime_salt(mmdd_indx{i}));
        reg_clim_do_b(i) = nanmean(regime_do(mmdd_indx{i}));
        reg_clim_no3_b(i) = nanmean(regime_no3(mmdd_indx{i}));
        reg_clim_temp_b(i) = nanmean(regime_temp(mmdd_indx{i}));
        reg_clim_salt_b(i) = nanmean(regime_salt(mmdd_indx{i}));
    end
end


%surface
%DO
clearvars xp_w_do pf_w_do
xp = 1:366;
xp_w_do = find(isnan(reg_clim_do)==0);
pf_w_do = polyfit(xp_w_do,reg_clim_do(xp_w_do),3);
yp_w_do_04 = polyval(pf_w_do,xp);
yp_w_do_04 = yp_w_do_04 + regm_do;

%NO3
clearvars xp_w_no3 pf_w_no3
xp = 1:366;
xp_w_no3 = find(isnan(reg_clim_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_clim_no3(xp_w_no3),3);
yp_w_no3_04 = polyval(pf_w_no3,xp);
yp_w_no3_04 = yp_w_no3_04 + regm_no3;

%salt
clearvars xp_w_salt pf_w_salt
xp = 1:366;
xp_w_salt = find(isnan(reg_clim_salt)==0);
pf_w_salt = polyfit(xp_w_salt,reg_clim_salt(xp_w_salt),3);
yp_w_salt_04 = polyval(pf_w_salt,xp);
yp_w_salt_04 = yp_w_salt_04 + regm_salt;

%temp
clearvars xp_w_temp pf_w_temp
xp = 1:366;
xp_w_temp = find(isnan(reg_clim_temp)==0);
pf_w_temp = polyfit(xp_w_temp,reg_clim_temp(xp_w_temp),3);
yp_w_temp_04 = polyval(pf_w_temp,xp);
yp_w_temp_04 = yp_w_temp_04 + regm_temp;


%bottom
%DO
clearvars xp_w_do_b pf_w_do_b
xp = 1:366;
xp_w_do_b = find(isnan(reg_clim_do_b)==0);
pf_w_do_b = polyfit(xp_w_do_b,reg_clim_do_b(xp_w_do_b),3);
yp_w_do_04_b = polyval(pf_w_do_b,xp);
yp_w_do_04_b = yp_w_do_04_b + regm_do_b;

%NO3
clearvars xp_w_no3_b pf_w_no3_b
xp = 1:366;
xp_w_no3_b = find(isnan(reg_clim_no3_b)==0);
pf_w_no3_b = polyfit(xp_w_no3_b,reg_clim_no3_b(xp_w_no3_b),3);
yp_w_no3_04_b = polyval(pf_w_no3_b,xp);
yp_w_no3_04_b = yp_w_no3_04_b + regm_no3_b;

%salt
clearvars xp_w_salt_b pf_w_salt_b
xp = 1:366;
xp_w_salt_b = find(isnan(reg_clim_salt_b)==0);
pf_w_salt_b = polyfit(xp_w_salt_b,reg_clim_salt_b(xp_w_salt_b),3);
yp_w_salt_04_b = polyval(pf_w_salt_b,xp);
yp_w_salt_04_b = yp_w_salt_04_b + regm_salt_b;

%temp
clearvars xp_w_temp_b pf_w_temp_b
xp = 1:366;
xp_w_temp_b = find(isnan(reg_clim_temp_b)==0);
pf_w_temp_b = polyfit(xp_w_temp_b,reg_clim_temp_b(xp_w_temp_b),3);
yp_w_temp_04_b = polyval(pf_w_temp_b,xp);
yp_w_temp_04_b = yp_w_temp_04_b + regm_temp_b;

return

%% DO
figure;
scatter(1:366,reg_clim_do + regm_do);
hold on
plot(1:366, yp_w_do_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('��������(ǥ��)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 18]) 
xlim([1 366])

%% salt
figure;
scatter(1:366,reg_clim_salt + regm_salt);
hold on
plot(1:366, yp_w_salt_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('salta (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial salt.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3 + regm_no3);
hold on
plot(1:366, yp_w_no3_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('��������(ǥ��)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% temp
figure;
scatter(1:366,reg_clim_temp + regm_temp);
hold on
plot(1:366, yp_w_temp_04,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('temp (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial temp.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])



%bottom
%% DO
figure;
scatter(1:366,reg_clim_do_b + regm_do_b);
hold on
plot(1:366, yp_w_do_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('DO (mg/L)','fontsize',13)
title('��������(����)-polynomial DO.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([2 18]) 
xlim([1 366])

%% salt
figure;
scatter(1:366,reg_clim_salt_b + regm_salt_b);
hold on
plot(1:366, yp_w_salt_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('salta (mg/m^3)','fontsize',13)
title('sumjin(songjung)-polynomial salt.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])
xlim([1 366])


%% no3
figure;
scatter(1:366,reg_clim_no3_b + regm_no3_b);
hold on
plot(1:366, yp_w_no3_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
title('��������(����)-polynomial no3.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% temp
figure;
scatter(1:366,reg_clim_temp_b + regm_temp_b);
hold on
plot(1:366, yp_w_temp_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('temp (mg/L)','fontsize',13)
title('sumjin(songjung)-polynomial temp.','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([0 inf])
xlim([1 366])

%% ALL
figure;
hold on;
% plot(1:366, yp_w_do_04_b,'--','color','b','linew',2)
% plot(1:366, yp_w_salt_04_b./1000,'--','color','g','linew',2)
plot(1:366, yp_w_no3_04,'--','color','b','linew',2)
plot(1:366, yp_w_no3_04_b,'--','color','r','linew',2)
xlabel('time(days)','fontsize',13)
ylabel('Concentration (ug/L)','fontsize',13)
title('��������(ǥ��&����)-polynomial.','fontsize',13)
grid on;
legend('no3-sur','no3-bot');
set(gca,'fontsize',13)
xlim([1 366])

%% plot raw
% make 1980~present
k=0
for i = 1980:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

clearvars t_tick
for i = 1:length(1980:2019)-1
    t_tick(i)=sum(t_tick_pre(1:i))+1;
end
t_tick=[1; t_tick';];
tx_tick = t_tick;

%plot-NO3
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= no3_sur;
   temp_b= no3_bot;
figure; hold on;
%  plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);   
   
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

%  plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 inf])
xlabel('time(year)','fontsize',13)
ylabel('no3 (ug/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:3:end)]);
set(gca,'xlim',[tx_tick(8) tx_tick(end)]);
set(gca,'xticklabel',1980:3:2019);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
% legend('surface','bottom')


%plot-DO
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= do_sur;
   temp_b= do_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);
%  plot(temp_sur,'*','color','c','linew',2);
%  plot(temp_bot,'*','color','m','linew',2);   
   
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
legend('surface','bottom')

%plot-temp
clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= temp_sur;
   temp_b= temp_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2);
 plot(temp_b,'*','color','r','linew',2);
%  plot(temp_sur,'*','color','c','linew',2);
%  plot(temp_bot,'*','color','m','linew',2);   
   
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

 plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([0 inf])
xlabel('time(year)','fontsize',13)
ylabel('temp (^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);
% set(gca,'Xminortick','on');
% set(gca,'Xminortickvalues',tx_tick(1:1:end));
legend('surface','bottom')


%plot- DO vs. temp
figure; hold on;
plot(temp_sur,do_sur,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
do_s_nonan=find(isnan(do_sur)==0);
do_b_nonan=find(isnan(do_bot)==0);

y_s_do_1 = polyfit(temp_sur(do_s_nonan), do_sur(do_s_nonan),1)
% re_do_sur = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_sur,'color','r');

yCalc = [5:0.01:30].* y_s_do_1(1) + y_s_do_1(2);
plot(5:0.01:30,yCalc,'-.','color','k','linew',2);
text(6,4.5,['y = ' num2str(y_s_do_1(1),'%.3f') 'x + ' num2str(y_s_do_1(2),'%.3f')],'fontsize',13);


%plot- DO vs. temp  bot
figure; hold on;
plot(temp_bot,do_bot,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
do_s_nonan=find(isnan(do_bot)==0);
do_b_nonan=find(isnan(do_bot)==0);
do_b_nonan(find(isnan(temp_bot(do_b_nonan)) == 1)) = [];

y_b_do_1 = polyfit(temp_bot(do_b_nonan), do_bot(do_b_nonan),1)
% re_do_bot = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_bot,'color','r');

yCalc_b = [0:0.01:30].* y_b_do_1(1) + y_b_do_1(2);
plot(0:0.01:30,yCalc_b,'-.','color','k','linew',2);
text(2,3,['y = ' num2str(y_b_do_1(1),'%.3f') 'x + ' num2str(y_b_do_1(2),'%.3f')],'fontsize',13);
plot(repmat([13],length(2:0.01:8),1),2:0.01:8,'-.','color','r','linew',2);


%plot-DO & recon.DO
% temp_s_nonan=find(isnan(temp_sur)==0);
% yCalc = temp_sur(temp_s_nonan) .* y_s_do_1(1) + y_s_do_1(2);

temp_b_nonan=find(isnan(temp_bot)==0);
yCalc_b = temp_bot(temp_b_nonan) .* y_b_do_1(1) + y_b_do_1(2);

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= do_sur;
   temp_b= do_bot;
   
figure;     
hold on;
%  plot(temp_s,'*','color','b','linew',2); 
% plot(temp_s_nonan, yCalc,'*','color','r','linew',2);
 plot(temp_b,'*','color','r','linew',2);
 plot(temp_b_nonan, yCalc_b,'*','color','b','linew',2);
    
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

%  plot(temp_s,'color','b');
 plot(temp_b,'color','r');  
 ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);

% legend('surface','surface_r_e')
legend('bottom','bottom_r_e')

set(gca,'xtick',[tx_tick(1:5:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:5:2019);
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [9, 7])

%plot- DO vs. DO  bot
figure; hold on;
plot(do_sur,do_bot,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('DO-sur (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
do_s_nonan=find(isnan(do_bot)==0);
do_b_nonan=find(isnan(do_bot)==0);
do_b_nonan(find(isnan(do_sur(do_b_nonan)) == 1)) = [];

y_b_do_2 = polyfit(do_sur(do_b_nonan), do_bot(do_b_nonan),1)
% re_do_bot = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_bot,'color','r');

yCalc_b_2 = [4:0.01:8].* y_b_do_2(1) + y_b_do_2(2);
plot(4:0.01:8,yCalc_b_2,'-.','color','k','linew',2);
text(4.5,7.3,['y = ' num2str(y_b_do_2(1),'%.3f') 'x + ' num2str(y_b_do_2(2),'%.3f')],'fontsize',13);



%plot-DO & recon.DO  2
% temp_s_nonan, yCalc
do_re_2 = yCalc .* y_b_do_2(1) + y_b_do_2(2);

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_b
   temp_b= do_bot;
figure; hold on;
 plot(temp_b,'*','color','r','linew',2);
 plot(temp_s_nonan, do_re_2,'*','color','b','linew',2);
    
t=1:length(no3_sur);
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

 plot(temp_b,'color','r');  
 ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);
legend('bottom','bottom_r_e')

set(gca,'xtick',[tx_tick(1:5:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:5:2019);

set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [9, 7])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%plot- DO vs. temps diff(surf-bot)
clearvars temp_diff*

below13=find(temp_bot < 13);
bt13=find(temp_bot >= 13);
temp_diff = temp_sur(bt13) - temp_bot(bt13);

% % make quasi_climate 1980~2019
% for i = 1:length(mmdd_indx)
%     if size(mmdd_indx{i},1) == 0     
%         reg_temp_diff{i} = {};
%     else
%         reg_temp_diff{i} = temp_diff(mmdd_indx{i})
%     end
% end
% 
% figure; hold on;
% for i = 1:366
%     plot(i,reg_temp_diff{i},'.','color','b');
% end
% xlabel('time(days)','fontsize',13)
% ylabel('temp diff(surf - bot) (^oC)','fontsize',13)
% title('temp diff(surf - bot)','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% ylim([-inf inf])
% xlim([1 366])
% 
% xx = repmat([1:366],40,1);
% clearvars xx_mm
% for i = 1:366
% if i == 1
%     xx_mm = xx(:,i);  
% elseif i~= 1
%     xx_mm = [xx_mm(:,1); xx(:,i);];
% end
% end
% 
% 
% size_c = [];
% for i = 1:366; size_c(i)=length(reg_temp_diff{i}); end
% bar(size_c)
% 
% clearvars temp_diff_mtx
% for i = 1:366
%     if length(reg_temp_diff{i}) == 40
%        temp_diff_mtx(:,i) = [reg_temp_diff{i}];  
%     elseif length(reg_temp_diff{i}) ~= 40
%        temp_diff_mtx(:,i) = [reg_temp_diff{i}; NaN(40 - length(reg_temp_diff{i}),1);];
%     end
% end
% 
% clearvars temp_diff_fix
% for i = 1:366
% if i == 1
%     temp_diff_fix = temp_diff_mtx(:,i);  
% elseif i~= 1
%     temp_diff_fix = [temp_diff_fix(:,1); temp_diff_mtx(:,i);];
% end
% end
% 
% diff_f_nan = find(isnan(temp_diff_fix)==1);
% xx_mm(diff_f_nan) = [];
% temp_diff_fix(diff_f_nan) = [];
% 
% plot(xx_mm,temp_diff_fix,'.','color','b');
% xlabel('time(days)','fontsize',13)
% ylabel('temp diff(surf - bot) (^oC)','fontsize',13)
% title('temp diff(surf - bot)','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% ylim([-inf inf])
% xlim([1 366])

%% bt13
figure; hold on;
plot(temp_diff, do_bot(bt13),'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
bt13_re=bt13;
bt13_re(find(isnan(do_bot(bt13)) == 1)) = []; %  do_bot(temp_diff_nonan) 's NaN removing
temp_diff_re =  temp_sur(bt13_re) - temp_bot(bt13_re);

y_b_do_diff = polyfit(temp_diff_re, do_bot(bt13_re),1)
% re_do_sur = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_sur,'color','r');

yCalc = [-2:0.01:16].* y_b_do_diff(1) + y_b_do_diff(2);
plot(-2:0.01:16,yCalc,'-.','color','k','linew',2);
text(6,6,['y = ' num2str(y_b_do_diff(1),'%.3f') 'x + ' num2str(y_b_do_diff(2),'%.3f')],'fontsize',13);

%%below 13
figure; hold on;
plot(temp_bot(below13), do_bot(below13),'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
below13_re=below13;
below13_re(find(isnan(do_bot(below13)) == 1)) = []; %  do_bot(temp_diff_nonan) 's NaN removing
temp_b_below = temp_bot(below13_re);

y_b_do_below = polyfit(temp_b_below, do_bot(below13_re),1)
% re_do_sur = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_sur,'color','r');

yCalc_below = [8:0.01:14].* y_b_do_below(1) + y_b_do_below(2);
plot(8:0.01:14,yCalc_below,'-.','color','k','linew',2);
text(9,4.5,['y = ' num2str(temp_b_below(1),'%.3f') 'x + ' num2str(temp_b_below(2),'%.3f')],'fontsize',13);


%% recon. diff T13
%plot-DO & recon.DO
clearvars temp_b yCalc_b_re yCal_con_int t t_in yCalc_below yCal_con_mer
% temp_b_nonan=find(isnan(temp_bot)==0);
yCalc_b_re = temp_diff .* y_b_do_diff(1) + y_b_do_diff(2);
yCalc_below = temp_b_below .* y_b_do_below(1) + y_b_do_below(2);

yCal_con_mer = [yCalc_b_re; yCalc_below;]; yCal_con_int = yCal_con_mer;
xx_con_mer = [bt13; below13_re;];

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_b= do_bot;
   
figure;     
hold on;
 plot(temp_b,'*','color','r','linew',2);
 plot(bt13, yCalc_b_re,'*','color','b','linew',2);
 plot(below13_re, yCalc_below, '*', 'color','c','linew',2);

 t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

[t t_in]=sort(xx_con_mer);
yCal_con_int = yCal_con_int(t_in);
yCal_con_int(isnan(yCal_con_int)) = interp1( t(~isnan(yCal_con_int)), yCal_con_int(~isnan(yCal_con_int)), t(isnan(yCal_con_int)) ); 

interp_nh4_b = temp_b;

 plot(temp_b,'color','r');
 plot(t,yCal_con_int,'color','g');  
 ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);

% legend('surface','surface_r_e')
legend('bottom_O_B_S','bottom_r_e_-_u_p_1_3','bottom_r_e_-_d_o_w_n_1_3')

set(gca,'xtick',[tx_tick(1:5:end)]);
set(gca,'xlim',[tx_tick(1) tx_tick(end)]);
set(gca,'xticklabel',1980:5:2019);
set(gcf, 'Units', 'Inches', 'Position', [1, 1, 9, 7], ...
    'PaperUnits', 'Inches', 'PaperSize', [9, 7])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% old %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%plot- DO vs. temps diff(surf-bot)
clearvars temp_diff*
temp_diff = temp_sur - temp_bot;

% make quasi_climate 1980~2019
for i = 1:length(mmdd_indx)
    if size(mmdd_indx{i},1) == 0     
        reg_temp_diff{i} = {};
    else
        reg_temp_diff{i} = temp_diff(mmdd_indx{i})
    end
end

figure; hold on;
for i = 1:366
    plot(i,reg_temp_diff{i},'.','color','b');
end
xlabel('time(days)','fontsize',13)
ylabel('temp diff(surf - bot) (^oC)','fontsize',13)
title('temp diff(surf - bot)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

xx = repmat([1:366],40,1);
clearvars xx_mm
for i = 1:366
if i == 1
    xx_mm = xx(:,i);  
elseif i~= 1
    xx_mm = [xx_mm(:,1); xx(:,i);];
end
end


size_c = [];
for i = 1:366; size_c(i)=length(reg_temp_diff{i}); end
bar(size_c)

clearvars temp_diff_mtx
for i = 1:366
    if length(reg_temp_diff{i}) == 40
       temp_diff_mtx(:,i) = [reg_temp_diff{i}];  
    elseif length(reg_temp_diff{i}) ~= 40
       temp_diff_mtx(:,i) = [reg_temp_diff{i}; NaN(40 - length(reg_temp_diff{i}),1);];
    end
end

clearvars temp_diff_fix
for i = 1:366
if i == 1
    temp_diff_fix = temp_diff_mtx(:,i);  
elseif i~= 1
    temp_diff_fix = [temp_diff_fix(:,1); temp_diff_mtx(:,i);];
end
end

diff_f_nan = find(isnan(temp_diff_fix)==1);
xx_mm(diff_f_nan) = [];
temp_diff_fix(diff_f_nan) = [];

plot(xx_mm,temp_diff_fix,'.','color','b');
xlabel('time(days)','fontsize',13)
ylabel('temp diff(surf - bot) (^oC)','fontsize',13)
title('temp diff(surf - bot)','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([-inf inf])
xlim([1 366])

figure; hold on;
plot(temp_diff, do_bot,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
temp_diff_nonan=find(isnan(temp_diff)==0);
temp_diff_nonan(find(isnan(do_bot(temp_diff_nonan)) == 1)) = []; %  do_bot(temp_diff_nonan) 's NaN removing


y_b_do_diff = polyfit(temp_diff(temp_diff_nonan), do_bot(temp_diff_nonan),1)
% re_do_sur = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_sur,'color','r');

yCalc = [-2:0.01:16].* y_b_do_diff(1) + y_b_do_diff(2);
plot(-2:0.01:16,yCalc,'-.','color','k','linew',2);
text(6,6.5,['y = ' num2str(y_s_do_1(1),'%.3f') 'x + ' num2str(y_s_do_1(2),'%.3f')],'fontsize',13);


%plot- DO vs. temp  bot
figure; hold on;
plot(temp_bot,do_bot,'.','color','b');
ylabel('DO (mg/L)','fontsize',13)
xlabel('Temp(^oC)','fontsize',13)
grid on; set(gca,'fontsize',13)

%regression
%slope y = b1*x
do_s_nonan=find(isnan(do_bot)==0);
do_b_nonan=find(isnan(do_bot)==0);
do_b_nonan(find(isnan(temp_bot(do_b_nonan)) == 1)) = [];

y_b_do_1 = polyfit(temp_bot(do_b_nonan), do_bot(do_b_nonan),1)
% re_do_bot = polyval(y_s_do_1, 5:0.01:30);
% plot(5:0.01:30,re_do_bot,'color','r');

yCalc_b = [0:0.01:30].* y_b_do_1(1) + y_b_do_1(2);
plot(0:0.01:30,yCalc_b,'-.','color','k','linew',2);
text(2,2,['y = ' num2str(y_b_do_1(1),'%.3f') 'x + ' num2str(y_b_do_1(2),'%.3f')],'fontsize',13);
