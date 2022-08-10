close all; clear; clc;

cd F:\ROMS\roms_tools\Run\
start

cd D:\장기생태\Dynamic\KOEM
% koem data
load koem_timeseires_monthly_3sig.mat
% location
load KOEM_st_info_(decimal_deg).mat
lat_koem = [lat_1; lat_2; lat_3;];
lon_koem = [lon_1; lon_2; lon_3;];

clearvars lat_1 lon_1 lat_2 lat_3 lon_2 lon_3
grd_file='grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
% get data from model result


% nearest point
for i = 1:length(lon_koem)
clearvars temp_2d near_temp
    temp_2d = sqrt((x - lon_koem(i)).^2 + (y - lat_koem(i)).^2);
    temp_2d = temp_2d .* (mask./mask);
    near_temp = find(nanmin(nanmin(temp_2d))==temp_2d);
    near_point(i) = near_temp(1);
    dist_2d(:,:,i) = temp_2d;
end

for i = 1:length(near_point)
    clearvars row col
    [row,col]=ind2sub(size(x),near_point(i)); %% 1d to 2d index
    near_point_2d(i,:) = [row, col]; 
end

% 65 points have some out of domain points
% find  it with this.
clearvars ind_out
ind_out=ones(1,65);
for i = 1:65
   if x(near_point_2d(i,1),near_point_2d(i,2)) > 128.15
    ind_out(i) = 0;
   end
   if y(near_point_2d(i,1),near_point_2d(i,2)) < 34.61
    ind_out(i) = 0;
   end
end


figure;
pcolor(x,y,mask); shading flat; hold on; 
for i = 1:65
    clearvars nanx
    nanx =  ind_out(i) ./ ind_out(i); %extract out of domain points;
    plot(x(near_point_2d(i,1),near_point_2d(i,2)).*nanx,y(near_point_2d(i,1),near_point_2d(i,2)).*nanx,'co');
end

% 49 month is jan.2001 (2,5,8,11 month was observed)
for i = 1:65
    if isnan(regime_no3(i,50)) == 1  
       ind_out(i) = 0;
    end
end

figure;
pcolor(x,y,mask); shading flat; hold on; 
for i = 1:65
    clearvars nanx
    nanx =  ind_out(i) ./ ind_out(i); %extract out of domain points;
    plot(x(near_point_2d(i,1),near_point_2d(i,2)).*nanx,y(near_point_2d(i,1),near_point_2d(i,2)).*nanx,'c+');
    text(x(near_point_2d(i,1),near_point_2d(i,2)).*nanx,y(near_point_2d(i,1),near_point_2d(i,2)).*nanx,['st-',num2str(i,'%02d')],'color','y');
end

% spatial region index
sp_gy = [4,22,28,29,30,32,33,34,35];
sp_s_gy = [38, 39, 42, 41, 44];
sp_gm = [6,9,10];% gamak bay
sp_e_gy = [48, 49];
sp_jj = [36, 37, 45, 46]; %jinju bay

% make eom_d
k=0
for i = 2001:2001
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end
t_tick_pre=sum(eom_d,2);

for i = 1:size(eom_d,1)
    for j = 1:size(eom_d,2)
        eom_d_each(i,j) = sum(eom_d(i,1:j));
    end
end
    
obs_tdx = 49:60;
for j=1:65 % num. of st.
for i=1:12 % month
    if i==1
       %surf
       obs_no3(j,1:eom_d_each(1,i)) = regime_no3(j,obs_tdx(i));
       obs_nh4(j,1:eom_d_each(1,i)) = regime_nh4(j,obs_tdx(i));
       obs_do(j,1:eom_d_each(1,i)) = regime_do(j,obs_tdx(i));
       obs_chl(j,1:eom_d_each(1,i)) = regime_chl(j,obs_tdx(i));
       obs_temp(j,1:eom_d_each(1,i)) = regime_temp(j,obs_tdx(i));
       obs_salt(j,1:eom_d_each(1,i)) = regime_salt(j,obs_tdx(i));
       % bot
       obs_no3_b(j,1:eom_d_each(1,i)) = regime_no3_b(j,obs_tdx(i));
       obs_nh4_b(j,1:eom_d_each(1,i)) = regime_nh4_b(j,obs_tdx(i));
       obs_do_b(j,1:eom_d_each(1,i)) = regime_do_b(j,obs_tdx(i));
       obs_chl_b(j,1:eom_d_each(1,i)) = regime_chl_b(j,obs_tdx(i));
       obs_temp_b(j,1:eom_d_each(1,i)) = regime_temp_b(j,obs_tdx(i));
       obs_salt_b(j,1:eom_d_each(1,i)) = regime_salt_b(j,obs_tdx(i));
    else
       obs_no3(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_no3(j,obs_tdx(i));
       obs_nh4(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_nh4(j,obs_tdx(i));
       obs_do(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_do(j,obs_tdx(i));
       obs_chl(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_chl(j,obs_tdx(i));
       obs_temp(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_temp(j,obs_tdx(i));
       obs_salt(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_salt(j,obs_tdx(i));
       % bot
       obs_no3_b(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_no3_b(j,obs_tdx(i));
       obs_nh4_b(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_nh4_b(j,obs_tdx(i));
       obs_do_b(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_do_b(j,obs_tdx(i));
       obs_chl_b(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_chl_b(j,obs_tdx(i));
       obs_temp_b(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_temp_b(j,obs_tdx(i));
       obs_salt_b(j,eom_d_each(1,i-1)+1:eom_d_each(1,i)) = regime_salt_b(j,obs_tdx(i));
    end
end
end

obs_tdx = 49:60;
for j=1:65 % num. of st.
for i=1:12 % month
       %surf
       clim_no3(j,i) = nanmean(regime_no3(j,i:12:end));
       clim_nh4(j,i) =  nanmean(regime_nh4(j,i:12:end));
       clim_do(j,i) =  nanmean(regime_do(j,i:12:end));
       clim_chl(j,i) =  nanmean(regime_chl(j,i:12:end));
       clim_temp(j,i) =  nanmean(regime_temp(j,i:12:end));
       clim_salt(j,i) =  nanmean(regime_salt(j,i:12:end));
       % bot
       clim_no3_b(j,i) =  nanmean(regime_no3_b(j,i:12:end));
       clim_nh4_b(j,i) =  nanmean(regime_nh4_b(j,i:12:end));
       clim_do_b(j,i) =  nanmean(regime_do_b(j,i:12:end));
       clim_chl_b(j,i) =  nanmean(regime_chl_b(j,i:12:end));
       clim_temp_b(j,i) =  nanmean(regime_temp_b(j,i:12:end));
       clim_salt_b(j,i) =  nanmean(regime_salt_b(j,i:12:end));
end
end





% compare plot KOEM vs. MODEL
n_2d_x = near_point_2d(:,1).*(ind_out'./ind_out'); n_2d_y = near_point_2d(:,2).*(ind_out'./ind_out');
n_2d_x(isnan(n_2d_x))=[]; n_2d_y(isnan(n_2d_y))=[];



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  spatial mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % spatial region index
sp_gy = [4,22,28,29,30,32,33,34,35];
sp_s_gy = [38, 39, 42, 41, 44];
sp_gm = [6,9,10];% gamak bay
sp_e_gy = [48, 49];
sp_jj = [36, 37, 45, 46]; %jinju bay

%% obs 
clearvars sp_std_*
for j = 1:length(sp_gy)
    clim_gy_no3(j,:)=squeeze(clim_no3(sp_gy(j),:));
    clim_gy_no3_b(j,:)=squeeze(clim_no3_b(sp_gy(j),:));
    clim_gy_nh4(j,:)=squeeze(clim_nh4(sp_gy(j),:));
    clim_gy_nh4_b(j,:)=squeeze(clim_nh4_b(sp_gy(j),:));
    clim_gy_chl(j,:)=squeeze(clim_chl(sp_gy(j),:));
    clim_gy_chl_b(j,:)=squeeze(clim_chl_b(sp_gy(j),:));
    clim_gy_temp(j,:)=squeeze(clim_temp(sp_gy(j),:));
    clim_gy_temp_b(j,:)=squeeze(clim_temp_b(sp_gy(j),:));
    clim_gy_salt(j,:)=squeeze(clim_salt(sp_gy(j),:));
    clim_gy_salt_b(j,:)=squeeze(clim_salt_b(sp_gy(j),:));
    clim_gy_do(j,:)=squeeze(clim_do(sp_gy(j),:));
    clim_gy_do_b(j,:)=squeeze(clim_do_b(sp_gy(j),:));  
end

for j = 1:length(sp_s_gy)
    clim_sgy_no3(j,:)=squeeze(clim_no3(sp_s_gy(j),:));
    clim_sgy_no3_b(j,:)=squeeze(clim_no3_b(sp_s_gy(j),:));
    clim_sgy_nh4(j,:)=squeeze(clim_nh4(sp_s_gy(j),:));
    clim_sgy_nh4_b(j,:)=squeeze(clim_nh4_b(sp_s_gy(j),:));
    clim_sgy_chl(j,:)=squeeze(clim_chl(sp_s_gy(j),:));
    clim_sgy_chl_b(j,:)=squeeze(clim_chl_b(sp_s_gy(j),:));
    clim_sgy_temp(j,:)=squeeze(clim_temp(sp_s_gy(j),:));
    clim_sgy_temp_b(j,:)=squeeze(clim_temp_b(sp_s_gy(j),:));
    clim_sgy_salt(j,:)=squeeze(clim_salt(sp_s_gy(j),:));
    clim_sgy_salt_b(j,:)=squeeze(clim_salt_b(sp_s_gy(j),:));
    clim_sgy_do(j,:)=squeeze(clim_do(sp_s_gy(j),:));
    clim_sgy_do_b(j,:)=squeeze(clim_do_b(sp_s_gy(j),:));  
end

for j = 1:length(sp_e_gy)
    clim_egy_no3(j,:)=squeeze(clim_no3(sp_e_gy(j),:));
    clim_egy_no3_b(j,:)=squeeze(clim_no3_b(sp_e_gy(j),:));
    clim_egy_nh4(j,:)=squeeze(clim_nh4(sp_e_gy(j),:));
    clim_egy_nh4_b(j,:)=squeeze(clim_nh4_b(sp_e_gy(j),:));
    clim_egy_chl(j,:)=squeeze(clim_chl(sp_e_gy(j),:));
    clim_egy_chl_b(j,:)=squeeze(clim_chl_b(sp_e_gy(j),:));
    clim_egy_temp(j,:)=squeeze(clim_temp(sp_e_gy(j),:));
    clim_egy_temp_b(j,:)=squeeze(clim_temp_b(sp_e_gy(j),:));
    clim_egy_salt(j,:)=squeeze(clim_salt(sp_e_gy(j),:));
    clim_egy_salt_b(j,:)=squeeze(clim_salt_b(sp_e_gy(j),:));
    clim_egy_do(j,:)=squeeze(clim_do(sp_e_gy(j),:));
    clim_egy_do_b(j,:)=squeeze(clim_do_b(sp_e_gy(j),:));  
end

for j = 1:length(sp_jj)
    clim_jj_no3(j,:)=squeeze(clim_no3(sp_jj(j),:));
    clim_jj_no3_b(j,:)=squeeze(clim_no3_b(sp_jj(j),:));
    clim_jj_nh4(j,:)=squeeze(clim_nh4(sp_jj(j),:));
    clim_jj_nh4_b(j,:)=squeeze(clim_nh4_b(sp_jj(j),:));
    clim_jj_chl(j,:)=squeeze(clim_chl(sp_jj(j),:));
    clim_jj_chl_b(j,:)=squeeze(clim_chl_b(sp_jj(j),:));
    clim_jj_temp(j,:)=squeeze(clim_temp(sp_jj(j),:));
    clim_jj_temp_b(j,:)=squeeze(clim_temp_b(sp_jj(j),:));
    clim_jj_salt(j,:)=squeeze(clim_salt(sp_jj(j),:));
    clim_jj_salt_b(j,:)=squeeze(clim_salt_b(sp_jj(j),:));
    clim_jj_do(j,:)=squeeze(clim_do(sp_jj(j),:));
    clim_jj_do_b(j,:)=squeeze(clim_do_b(sp_jj(j),:));  
end


for i = 1:12
    %std
    clim_std_gy_no3(i)=nanstd(squeeze(clim_gy_no3(:,i)));
    clim_std_gy_no3_b(i)=nanstd(squeeze(clim_gy_no3_b(:,i)));
    clim_std_gy_nh4(i)=nanstd(squeeze(clim_gy_nh4(:,i)));
    clim_std_gy_nh4_b(i)=nanstd(squeeze(clim_gy_nh4_b(:,i)));
    clim_std_gy_chl(i)=nanstd(squeeze(clim_gy_chl(:,i)));
    clim_std_gy_chl_b(i)=nanstd(squeeze(clim_gy_chl_b(:,i)));
    clim_std_gy_temp(i)=nanstd(squeeze(clim_gy_temp(:,i)));
    clim_std_gy_temp_b(i)=nanstd(squeeze(clim_gy_temp_b(:,i)));
    clim_std_gy_salt(i)=nanstd(squeeze(clim_gy_salt(:,i)));
    clim_std_gy_salt_b(i)=nanstd(squeeze(clim_gy_salt_b(:,i)));
    clim_std_gy_do(i)=nanstd(squeeze(clim_gy_do(:,i)));
    clim_std_gy_do_b(i)=nanstd(squeeze(clim_gy_do_b(:,i)));

    clim_std_sgy_no3(i)=nanstd(squeeze(clim_sgy_no3(:,i)));
    clim_std_sgy_no3_b(i)=nanstd(squeeze(clim_sgy_no3_b(:,i)));
    clim_std_sgy_nh4(i)=nanstd(squeeze(clim_sgy_nh4(:,i)));
    clim_std_sgy_nh4_b(i)=nanstd(squeeze(clim_sgy_nh4_b(:,i)));
    clim_std_sgy_chl(i)=nanstd(squeeze(clim_sgy_chl(:,i)));
    clim_std_sgy_chl_b(i)=nanstd(squeeze(clim_sgy_chl_b(:,i)));
    clim_std_sgy_temp(i)=nanstd(squeeze(clim_sgy_temp(:,i)));
    clim_std_sgy_temp_b(i)=nanstd(squeeze(clim_sgy_temp_b(:,i)));
    clim_std_sgy_salt(i)=nanstd(squeeze(clim_sgy_salt(:,i)));
    clim_std_sgy_salt_b(i)=nanstd(squeeze(clim_sgy_salt_b(:,i)));
    clim_std_sgy_do(i)=nanstd(squeeze(clim_sgy_do(:,i)));
    clim_std_sgy_do_b(i)=nanstd(squeeze(clim_sgy_do_b(:,i)));

    clim_std_egy_no3(i)=nanstd(squeeze(clim_egy_no3(:,i)));
    clim_std_egy_no3_b(i)=nanstd(squeeze(clim_egy_no3_b(:,i)));
    clim_std_egy_nh4(i)=nanstd(squeeze(clim_egy_nh4(:,i)));
    clim_std_egy_nh4_b(i)=nanstd(squeeze(clim_egy_nh4_b(:,i)));
    clim_std_egy_chl(i)=nanstd(squeeze(clim_egy_chl(:,i)));
    clim_std_egy_chl_b(i)=nanstd(squeeze(clim_egy_chl_b(:,i)));
    clim_std_egy_temp(i)=nanstd(squeeze(clim_egy_temp(:,i)));
    clim_std_egy_temp_b(i)=nanstd(squeeze(clim_egy_temp_b(:,i)));
    clim_std_egy_salt(i)=nanstd(squeeze(clim_egy_salt(:,i)));
    clim_std_egy_salt_b(i)=nanstd(squeeze(clim_egy_salt_b(:,i)));
    clim_std_egy_do(i)=nanstd(squeeze(clim_egy_do(:,i)));
    clim_std_egy_do_b(i)=nanstd(squeeze(clim_egy_do_b(:,i)));

    clim_std_jj_no3(i)=nanstd(squeeze(clim_jj_no3(:,i)));
    clim_std_jj_no3_b(i)=nanstd(squeeze(clim_jj_no3_b(:,i)));
    clim_std_jj_nh4(i)=nanstd(squeeze(clim_jj_nh4(:,i)));
    clim_std_jj_nh4_b(i)=nanstd(squeeze(clim_jj_nh4_b(:,i)));
    clim_std_jj_chl(i)=nanstd(squeeze(clim_jj_chl(:,i)));
    clim_std_jj_chl_b(i)=nanstd(squeeze(clim_jj_chl_b(:,i)));
    clim_std_jj_temp(i)=nanstd(squeeze(clim_jj_temp(:,i)));
    clim_std_jj_temp_b(i)=nanstd(squeeze(clim_jj_temp_b(:,i)));
    clim_std_jj_salt(i)=nanstd(squeeze(clim_jj_salt(:,i)));
    clim_std_jj_salt_b(i)=nanstd(squeeze(clim_jj_salt_b(:,i)));
    clim_std_jj_do(i)=nanstd(squeeze(clim_jj_do(:,i)));
    clim_std_jj_do_b(i)=nanstd(squeeze(clim_jj_do_b(:,i)));
    
    %mean
    obm_gy_no3(i)=nanmean(squeeze(clim_gy_no3(:,i)));
    obm_gy_no3_b(i)=nanmean(squeeze(clim_gy_no3_b(:,i)));
    obm_gy_nh4(i)=nanmean(squeeze(clim_gy_nh4(:,i)));
    obm_gy_nh4_b(i)=nanmean(squeeze(clim_gy_nh4_b(:,i)));
    obm_gy_chl(i)=nanmean(squeeze(clim_gy_chl(:,i)));
    obm_gy_chl_b(i)=nanmean(squeeze(clim_gy_chl_b(:,i)));
    obm_gy_temp(i)=nanmean(squeeze(clim_gy_temp(:,i)));
    obm_gy_temp_b(i)=nanmean(squeeze(clim_gy_temp_b(:,i)));
    obm_gy_salt(i)=nanmean(squeeze(clim_gy_salt(:,i)));
    obm_gy_salt_b(i)=nanmean(squeeze(clim_gy_salt_b(:,i)));
    obm_gy_do(i)=nanmean(squeeze(clim_gy_do(:,i)));
    obm_gy_do_b(i)=nanmean(squeeze(clim_gy_do_b(:,i)));

    obm_sgy_no3(i)=nanmean(squeeze(clim_sgy_no3(:,i)));
    obm_sgy_no3_b(i)=nanmean(squeeze(clim_sgy_no3_b(:,i)));
    obm_sgy_nh4(i)=nanmean(squeeze(clim_sgy_nh4(:,i)));
    obm_sgy_nh4_b(i)=nanmean(squeeze(clim_sgy_nh4_b(:,i)));
    obm_sgy_chl(i)=nanmean(squeeze(clim_sgy_chl(:,i)));
    obm_sgy_chl_b(i)=nanmean(squeeze(clim_sgy_chl_b(:,i)));
    obm_sgy_temp(i)=nanmean(squeeze(clim_sgy_temp(:,i)));
    obm_sgy_temp_b(i)=nanmean(squeeze(clim_sgy_temp_b(:,i)));
    obm_sgy_salt(i)=nanmean(squeeze(clim_sgy_salt(:,i)));
    obm_sgy_salt_b(i)=nanmean(squeeze(clim_sgy_salt_b(:,i)));
    obm_sgy_do(i)=nanmean(squeeze(clim_sgy_do(:,i)));
    obm_sgy_do_b(i)=nanmean(squeeze(clim_sgy_do_b(:,i)));

    obm_egy_no3(i)=nanmean(squeeze(clim_egy_no3(:,i)));
    obm_egy_no3_b(i)=nanmean(squeeze(clim_egy_no3_b(:,i)));
    obm_egy_nh4(i)=nanmean(squeeze(clim_egy_nh4(:,i)));
    obm_egy_nh4_b(i)=nanmean(squeeze(clim_egy_nh4_b(:,i)));
    obm_egy_chl(i)=nanmean(squeeze(clim_egy_chl(:,i)));
    obm_egy_chl_b(i)=nanmean(squeeze(clim_egy_chl_b(:,i)));
    obm_egy_temp(i)=nanmean(squeeze(clim_egy_temp(:,i)));
    obm_egy_temp_b(i)=nanmean(squeeze(clim_egy_temp_b(:,i)));
    obm_egy_salt(i)=nanmean(squeeze(clim_egy_salt(:,i)));
    obm_egy_salt_b(i)=nanmean(squeeze(clim_egy_salt_b(:,i)));
    obm_egy_do(i)=nanmean(squeeze(clim_egy_do(:,i)));
    obm_egy_do_b(i)=nanmean(squeeze(clim_egy_do_b(:,i)));

    obm_jj_no3(i)=nanmean(squeeze(clim_jj_no3(:,i)));
    obm_jj_no3_b(i)=nanmean(squeeze(clim_jj_no3_b(:,i)));
    obm_jj_nh4(i)=nanmean(squeeze(clim_jj_nh4(:,i)));
    obm_jj_nh4_b(i)=nanmean(squeeze(clim_jj_nh4_b(:,i)));
    obm_jj_chl(i)=nanmean(squeeze(clim_jj_chl(:,i)));
    obm_jj_chl_b(i)=nanmean(squeeze(clim_jj_chl_b(:,i)));
    obm_jj_temp(i)=nanmean(squeeze(clim_jj_temp(:,i)));
    obm_jj_temp_b(i)=nanmean(squeeze(clim_jj_temp_b(:,i)));
    obm_jj_salt(i)=nanmean(squeeze(clim_jj_salt(:,i)));
    obm_jj_salt_b(i)=nanmean(squeeze(clim_jj_salt_b(:,i)));
    obm_jj_do(i)=nanmean(squeeze(clim_jj_do(:,i)));
    obm_jj_do_b(i)=nanmean(squeeze(clim_jj_do_b(:,i)));

end

return

%% GY plot

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
% gy %
errorbar(1:12, obm_gy_no3./14, clim_std_gy_no3./14,'color','r','linew',2); xlim([1 12]);
plot(obm_gy_no3./14,'kh','linew',2); 
% south boundary %
errorbar(1:12, obm_sgy_no3./14, clim_std_sgy_no3./14,'color','b','linew',2);      
plot(obm_sgy_no3./14,'rh','linew',2); 
% east boundary %
errorbar(1:12, obm_egy_no3./14, clim_std_egy_no3./14,'color','g','linew',2);      
plot(obm_egy_no3./14,'bh','linew',2); 

title(['GY KOEM NO3 spatial mean climate']);
xlabel('time(month)','fontsize',13)
ylabel('NO3 (umol/m^3)','fontsize',13)
grid on
set(gca,'fontsize',13)
print(fig,strcat('KOEM_spatial_mean_climate_no3'),'-dpng')


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
% gy %
errorbar(1:12, obm_gy_nh4./14, clim_std_gy_nh4./14,'color','r','linew',2); xlim([1 12]);
plot(obm_gy_nh4./14,'kh','linew',2); 
% south boundary %
errorbar(1:12, obm_sgy_nh4./14, clim_std_sgy_nh4./14,'color','b','linew',2);      
plot(obm_sgy_nh4./14,'rh','linew',2); 
% east boundary %
errorbar(1:12, obm_egy_nh4./14, clim_std_egy_nh4./14,'color','g','linew',2);      
plot(obm_egy_nh4./14,'bh','linew',2); 
% jinju bay %
errorbar(1:12, obm_jj_nh4./14, clim_std_jj_nh4./14,'color','m','linew',2);      
plot(obm_jj_nh4./14,'ch','linew',2); 

title(['GY KOEM NH4 spatial mean climate']);
xlabel('time(month)','fontsize',13)
ylabel('NO3 (umol/m^3)','fontsize',13)
grid on
set(gca,'fontsize',13)
print(fig,strcat('KOEM_spatial_mean_climate_nh4'),'-dpng')

%nh4
obm_gy_nh4(isnan(obm_gy_nh4)) = [];
obm_egy_nh4(isnan(obm_egy_nh4)) = [];
obm_sgy_nh4(isnan(obm_sgy_nh4)) = [];
obm_jj_nh4(isnan(obm_jj_nh4)) = [];

clim_std_gy_nh4(isnan(clim_std_gy_nh4)) = [];
clim_std_egy_nh4(isnan(clim_std_egy_nh4)) = [];
clim_std_sgy_nh4(isnan(clim_std_sgy_nh4)) = [];
clim_std_jj_nh4(isnan(clim_std_jj_nh4)) = [];

y=[obm_gy_nh4; obm_sgy_nh4; obm_egy_nh4; obm_jj_nh4;]
err=[clim_std_gy_nh4; clim_std_sgy_nh4; clim_std_egy_nh4; clim_std_jj_nh4;];

fig=figure; 
errorbar(y'./14,err'./14,'linew',2);
le=legend('GY','South','East','Jinju');
le.Location = 'northwest'
le.Orientation = 'horizontal'
title(['GY KOEM NH4 spatial mean climate']);
ylabel('NO3 (umol/m^3)','fontsize',13)
xlim([0.75 4.25])
xticks([1 2 3 4] );
xticklabels({'Feb','May','Aug','Nov'});
grid on
set(gca,'fontsize',13)
print(fig,strcat('KOEM_spatial_mean_climate_nh4'),'-dpng')

% no3
obm_gy_no3(isnan(obm_gy_no3)) = [];
obm_egy_no3(isnan(obm_egy_no3)) = [];
obm_sgy_no3(isnan(obm_sgy_no3)) = [];
obm_jj_no3(isnan(obm_jj_no3)) = [];

clim_std_gy_no3(isnan(clim_std_gy_no3)) = [];
clim_std_egy_no3(isnan(clim_std_egy_no3)) = [];
clim_std_sgy_no3(isnan(clim_std_sgy_no3)) = [];
clim_std_jj_no3(isnan(clim_std_jj_no3)) = [];

fig = figure; hold on;
y=[obm_gy_no3; obm_sgy_no3; obm_egy_no3; obm_jj_no3;]
err=[clim_std_gy_no3; clim_std_sgy_no3; clim_std_egy_no3; clim_std_jj_no3;];

fig=figure; 
errorbar(y'./14,err'./14,'linew',2);
le=legend('GY','South','East','Jinju');
le.Location = 'northwest'
le.Orientation = 'horizontal'
title(['GY KOEM no3 spatial mean climate']);
ylabel('NO3 (umol/m^3)','fontsize',13)
xlim([0.75 4.25])
xticks([1 2 3 4] );
xticklabels({'Feb','May','Aug','Nov'});
grid on
set(gca,'fontsize',13)
print(fig,strcat('KOEM_spatial_mean_climate_no3'),'-dpng')

%% KODC
% save('KODC_data_monthly_v5_40016.mat','*_bot','*_sur','ref_yymm', 'date_ymd');
% save('KODC_data_monthly_v5_20501.mat','*_bot','*_sur','ref_yymm', 'date_ymd');

kodc_400=load('KODC_data_monthly_v5_40016.mat');
kodc_205=load('KODC_data_monthly_v5_20501.mat');


kodc_no3_bot = [kodc_400.no3_bot; kodc_205.no3_bot;];
kodc_no3_surf = [kodc_400.no3_sur; kodc_205.no3_sur;];
kodc_205.no3_bot

for j = 1:2
for i=1:12 % month
       %surf
       kodc_clim_no3_bot(j,i) = nanmean(kodc_no3_bot(j,i:12:end));
       kodc_clim_no3_surf(j,i) =  nanmean(kodc_no3_surf(j,i:12:end));
end
end

for i = 1:12
    %std
    clim_std_kodc_no3(i)=nanstd(squeeze(kodc_clim_no3_surf(:,i)));
    clim_std_kodc_no3_b(i)=nanstd(squeeze(kodc_clim_no3_bot(:,i)));    
    %mean
    obm_kodc_no3(i)=nanmean(squeeze(kodc_clim_no3_surf(:,i)));
    obm_kodc_no3_b(i)=nanmean(squeeze(kodc_clim_no3_bot(:,i)));  
end

inx = [2,4,6,8,10,12];
clim_std_kodc_no3_in = clim_std_kodc_no3(inx);
clim_std_kodc_no3_b_in = clim_std_kodc_no3_b(inx);

obm_kodc_no3_in = obm_kodc_no3(inx);
obm_kodc_no3_b_in = obm_kodc_no3_b(inx);

fig_x =[2,5,8,11]';
for i = 2:4
fig_x(:,i) = [2,5,8,11]';
end

fig=figure; hold on;
errorbar(fig_x,y'./14,err'./14,'linew',2);
title(['GY KOEM no3 spatial mean climate']);
xlabel('month','fontsize',13)
ylabel('NO3 (umol/m^3)','fontsize',13)
xlim([1 12])
xticks([2 4 5 6 8 10 11 12] );
% xticklabels({'Feb','May','Aug','Nov'});
grid on
errorbar(inx,obm_kodc_no3_in',obm_kodc_no3_in','linew',2);
le=legend('GY','South','East','Jinju','KODC');
le.Location = 'northwest'
le.Orientation = 'horizontal'
set(gca,'fontsize',13)
print(fig,strcat('KOEM&KODC_spatial_mean_climate_no3'),'-dpng')

save('plot_KOEM_st_coastal_vs_remote_nitrate_data.mat','-v7.3');


