% spatial region index
sp_gy = [4,22,28,29,30,32,33,34,35];
sp_s_gy = [38, 39, 42, 41, 44];
sp_gm = [6,9,10];% gamak bay
sp_e_gy = [48, 49];
sp_jj = [36, 37, 45, 46]; %jinju bay

squeeze(nanmean(nanmean(mod_no3(near_point_2d(sp_gy,1), near_point_2d(sp_gy,2),20,:),1),2));
clearvars sp_std_*
for j = 1:length(sp_gy)
    sp_gy_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_do(j,:)=squeeze(mod_do(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));  
end

for j = 1:length(sp_s_gy)
    sp_sgy_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_do(j,:)=squeeze(mod_do(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));  
end

for j = 1:length(sp_e_gy)
    sp_egy_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_do(j,:)=squeeze(mod_do(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));  
end

for j = 1:length(sp_jj)
    sp_jj_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_do(j,:)=squeeze(mod_do(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));  
end


for i = 1:366
    %std
    std_gy_no3(i)=nanstd(squeeze(sp_gy_no3(:,i)));
    std_gy_no3_b(i)=nanstd(squeeze(sp_gy_no3_b(:,i)));
    std_gy_nh4(i)=nanstd(squeeze(sp_gy_nh4(:,i)));
    std_gy_nh4_b(i)=nanstd(squeeze(sp_gy_nh4_b(:,i)));
    std_gy_chl(i)=nanstd(squeeze(sp_gy_chl(:,i)));
    std_gy_chl_b(i)=nanstd(squeeze(sp_gy_chl_b(:,i)));
    std_gy_temp(i)=nanstd(squeeze(sp_gy_temp(:,i)));
    std_gy_temp_b(i)=nanstd(squeeze(sp_gy_temp_b(:,i)));
    std_gy_salt(i)=nanstd(squeeze(sp_gy_salt(:,i)));
    std_gy_salt_b(i)=nanstd(squeeze(sp_gy_salt_b(:,i)));
    std_gy_do(i)=nanstd(squeeze(sp_gy_do(:,i)));
    std_gy_do_b(i)=nanstd(squeeze(sp_gy_do_b(:,i)));

    std_sgy_no3(i)=nanstd(squeeze(sp_sgy_no3(:,i)));
    std_sgy_no3_b(i)=nanstd(squeeze(sp_sgy_no3_b(:,i)));
    std_sgy_nh4(i)=nanstd(squeeze(sp_sgy_nh4(:,i)));
    std_sgy_nh4_b(i)=nanstd(squeeze(sp_sgy_nh4_b(:,i)));
    std_sgy_chl(i)=nanstd(squeeze(sp_sgy_chl(:,i)));
    std_sgy_chl_b(i)=nanstd(squeeze(sp_sgy_chl_b(:,i)));
    std_sgy_temp(i)=nanstd(squeeze(sp_sgy_temp(:,i)));
    std_sgy_temp_b(i)=nanstd(squeeze(sp_sgy_temp_b(:,i)));
    std_sgy_salt(i)=nanstd(squeeze(sp_sgy_salt(:,i)));
    std_sgy_salt_b(i)=nanstd(squeeze(sp_sgy_salt_b(:,i)));
    std_sgy_do(i)=nanstd(squeeze(sp_sgy_do(:,i)));
    std_sgy_do_b(i)=nanstd(squeeze(sp_sgy_do_b(:,i)));

    std_egy_no3(i)=nanstd(squeeze(sp_egy_no3(:,i)));
    std_egy_no3_b(i)=nanstd(squeeze(sp_egy_no3_b(:,i)));
    std_egy_nh4(i)=nanstd(squeeze(sp_egy_nh4(:,i)));
    std_egy_nh4_b(i)=nanstd(squeeze(sp_egy_nh4_b(:,i)));
    std_egy_chl(i)=nanstd(squeeze(sp_egy_chl(:,i)));
    std_egy_chl_b(i)=nanstd(squeeze(sp_egy_chl_b(:,i)));
    std_egy_temp(i)=nanstd(squeeze(sp_egy_temp(:,i)));
    std_egy_temp_b(i)=nanstd(squeeze(sp_egy_temp_b(:,i)));
    std_egy_salt(i)=nanstd(squeeze(sp_egy_salt(:,i)));
    std_egy_salt_b(i)=nanstd(squeeze(sp_egy_salt_b(:,i)));
    std_egy_do(i)=nanstd(squeeze(sp_egy_do(:,i)));
    std_egy_do_b(i)=nanstd(squeeze(sp_egy_do_b(:,i)));

    std_jj_no3(i)=nanstd(squeeze(sp_jj_no3(:,i)));
    std_jj_no3_b(i)=nanstd(squeeze(sp_jj_no3_b(:,i)));
    std_jj_nh4(i)=nanstd(squeeze(sp_jj_nh4(:,i)));
    std_jj_nh4_b(i)=nanstd(squeeze(sp_jj_nh4_b(:,i)));
    std_jj_chl(i)=nanstd(squeeze(sp_jj_chl(:,i)));
    std_jj_chl_b(i)=nanstd(squeeze(sp_jj_chl_b(:,i)));
    std_jj_temp(i)=nanstd(squeeze(sp_jj_temp(:,i)));
    std_jj_temp_b(i)=nanstd(squeeze(sp_jj_temp_b(:,i)));
    std_jj_salt(i)=nanstd(squeeze(sp_jj_salt(:,i)));
    std_jj_salt_b(i)=nanstd(squeeze(sp_jj_salt_b(:,i)));
    std_jj_do(i)=nanstd(squeeze(sp_jj_do(:,i)));
    std_jj_do_b(i)=nanstd(squeeze(sp_jj_do_b(:,i)));
    
    %mean
    gy_no3(i)=nanmean(squeeze(sp_gy_no3(:,i)));
    gy_no3_b(i)=nanmean(squeeze(sp_gy_no3_b(:,i)));
    gy_nh4(i)=nanmean(squeeze(sp_gy_nh4(:,i)));
    gy_nh4_b(i)=nanmean(squeeze(sp_gy_nh4_b(:,i)));
    gy_chl(i)=nanmean(squeeze(sp_gy_chl(:,i)));
    gy_chl_b(i)=nanmean(squeeze(sp_gy_chl_b(:,i)));
    gy_temp(i)=nanmean(squeeze(sp_gy_temp(:,i)));
    gy_temp_b(i)=nanmean(squeeze(sp_gy_temp_b(:,i)));
    gy_salt(i)=nanmean(squeeze(sp_gy_salt(:,i)));
    gy_salt_b(i)=nanmean(squeeze(sp_gy_salt_b(:,i)));
    gy_do(i)=nanmean(squeeze(sp_gy_do(:,i)));
    gy_do_b(i)=nanmean(squeeze(sp_gy_do_b(:,i)));

    sgy_no3(i)=nanmean(squeeze(sp_sgy_no3(:,i)));
    sgy_no3_b(i)=nanmean(squeeze(sp_sgy_no3_b(:,i)));
    sgy_nh4(i)=nanmean(squeeze(sp_sgy_nh4(:,i)));
    sgy_nh4_b(i)=nanmean(squeeze(sp_sgy_nh4_b(:,i)));
    sgy_chl(i)=nanmean(squeeze(sp_sgy_chl(:,i)));
    sgy_chl_b(i)=nanmean(squeeze(sp_sgy_chl_b(:,i)));
    sgy_temp(i)=nanmean(squeeze(sp_sgy_temp(:,i)));
    sgy_temp_b(i)=nanmean(squeeze(sp_sgy_temp_b(:,i)));
    sgy_salt(i)=nanmean(squeeze(sp_sgy_salt(:,i)));
    sgy_salt_b(i)=nanmean(squeeze(sp_sgy_salt_b(:,i)));
    sgy_do(i)=nanmean(squeeze(sp_sgy_do(:,i)));
    sgy_do_b(i)=nanmean(squeeze(sp_sgy_do_b(:,i)));

    egy_no3(i)=nanmean(squeeze(sp_egy_no3(:,i)));
    egy_no3_b(i)=nanmean(squeeze(sp_egy_no3_b(:,i)));
    egy_nh4(i)=nanmean(squeeze(sp_egy_nh4(:,i)));
    egy_nh4_b(i)=nanmean(squeeze(sp_egy_nh4_b(:,i)));
    egy_chl(i)=nanmean(squeeze(sp_egy_chl(:,i)));
    egy_chl_b(i)=nanmean(squeeze(sp_egy_chl_b(:,i)));
    egy_temp(i)=nanmean(squeeze(sp_egy_temp(:,i)));
    egy_temp_b(i)=nanmean(squeeze(sp_egy_temp_b(:,i)));
    egy_salt(i)=nanmean(squeeze(sp_egy_salt(:,i)));
    egy_salt_b(i)=nanmean(squeeze(sp_egy_salt_b(:,i)));
    egy_do(i)=nanmean(squeeze(sp_egy_do(:,i)));
    egy_do_b(i)=nanmean(squeeze(sp_egy_do_b(:,i)));

    jj_no3(i)=nanmean(squeeze(sp_jj_no3(:,i)));
    jj_no3_b(i)=nanmean(squeeze(sp_jj_no3_b(:,i)));
    jj_nh4(i)=nanmean(squeeze(sp_jj_nh4(:,i)));
    jj_nh4_b(i)=nanmean(squeeze(sp_jj_nh4_b(:,i)));
    jj_chl(i)=nanmean(squeeze(sp_jj_chl(:,i)));
    jj_chl_b(i)=nanmean(squeeze(sp_jj_chl_b(:,i)));
    jj_temp(i)=nanmean(squeeze(sp_jj_temp(:,i)));
    jj_temp_b(i)=nanmean(squeeze(sp_jj_temp_b(:,i)));
    jj_salt(i)=nanmean(squeeze(sp_jj_salt(:,i)));
    jj_salt_b(i)=nanmean(squeeze(sp_jj_salt_b(:,i)));
    jj_do(i)=nanmean(squeeze(sp_jj_do(:,i)));
    jj_do_b(i)=nanmean(squeeze(sp_jj_do_b(:,i)));

end


%% obs 
clearvars sp_std_*
for j = 1:length(sp_gy)
    obs_gy_no3(j,:)=squeeze(obs_no3(sp_gy(j),:));
    obs_gy_no3_b(j,:)=squeeze(obs_no3_b(sp_gy(j),:));
    obs_gy_nh4(j,:)=squeeze(obs_nh4(sp_gy(j),:));
    obs_gy_nh4_b(j,:)=squeeze(obs_nh4_b(sp_gy(j),:));
    obs_gy_chl(j,:)=squeeze(obs_chl(sp_gy(j),:));
    obs_gy_chl_b(j,:)=squeeze(obs_chl_b(sp_gy(j),:));
    obs_gy_temp(j,:)=squeeze(obs_temp(sp_gy(j),:));
    obs_gy_temp_b(j,:)=squeeze(obs_temp_b(sp_gy(j),:));
    obs_gy_salt(j,:)=squeeze(obs_salt(sp_gy(j),:));
    obs_gy_salt_b(j,:)=squeeze(obs_salt_b(sp_gy(j),:));
    obs_gy_do(j,:)=squeeze(obs_do(sp_gy(j),:));
    obs_gy_do_b(j,:)=squeeze(obs_do_b(sp_gy(j),:));  
end

for j = 1:length(sp_s_gy)
    obs_sgy_no3(j,:)=squeeze(obs_no3(sp_s_gy(j),:));
    obs_sgy_no3_b(j,:)=squeeze(obs_no3_b(sp_s_gy(j),:));
    obs_sgy_nh4(j,:)=squeeze(obs_nh4(sp_s_gy(j),:));
    obs_sgy_nh4_b(j,:)=squeeze(obs_nh4_b(sp_s_gy(j),:));
    obs_sgy_chl(j,:)=squeeze(obs_chl(sp_s_gy(j),:));
    obs_sgy_chl_b(j,:)=squeeze(obs_chl_b(sp_s_gy(j),:));
    obs_sgy_temp(j,:)=squeeze(obs_temp(sp_s_gy(j),:));
    obs_sgy_temp_b(j,:)=squeeze(obs_temp_b(sp_s_gy(j),:));
    obs_sgy_salt(j,:)=squeeze(obs_salt(sp_s_gy(j),:));
    obs_sgy_salt_b(j,:)=squeeze(obs_salt_b(sp_s_gy(j),:));
    obs_sgy_do(j,:)=squeeze(obs_do(sp_s_gy(j),:));
    obs_sgy_do_b(j,:)=squeeze(obs_do_b(sp_s_gy(j),:));  
end

for j = 1:length(sp_e_gy)
    obs_egy_no3(j,:)=squeeze(obs_no3(sp_e_gy(j),:));
    obs_egy_no3_b(j,:)=squeeze(obs_no3_b(sp_e_gy(j),:));
    obs_egy_nh4(j,:)=squeeze(obs_nh4(sp_e_gy(j),:));
    obs_egy_nh4_b(j,:)=squeeze(obs_nh4_b(sp_e_gy(j),:));
    obs_egy_chl(j,:)=squeeze(obs_chl(sp_e_gy(j),:));
    obs_egy_chl_b(j,:)=squeeze(obs_chl_b(sp_e_gy(j),:));
    obs_egy_temp(j,:)=squeeze(obs_temp(sp_e_gy(j),:));
    obs_egy_temp_b(j,:)=squeeze(obs_temp_b(sp_e_gy(j),:));
    obs_egy_salt(j,:)=squeeze(obs_salt(sp_e_gy(j),:));
    obs_egy_salt_b(j,:)=squeeze(obs_salt_b(sp_e_gy(j),:));
    obs_egy_do(j,:)=squeeze(obs_do(sp_e_gy(j),:));
    obs_egy_do_b(j,:)=squeeze(obs_do_b(sp_e_gy(j),:));  
end

for j = 1:length(sp_jj)
    obs_jj_no3(j,:)=squeeze(obs_no3(sp_jj(j),:));
    obs_jj_no3_b(j,:)=squeeze(obs_no3_b(sp_jj(j),:));
    obs_jj_nh4(j,:)=squeeze(obs_nh4(sp_jj(j),:));
    obs_jj_nh4_b(j,:)=squeeze(obs_nh4_b(sp_jj(j),:));
    obs_jj_chl(j,:)=squeeze(obs_chl(sp_jj(j),:));
    obs_jj_chl_b(j,:)=squeeze(obs_chl_b(sp_jj(j),:));
    obs_jj_temp(j,:)=squeeze(obs_temp(sp_jj(j),:));
    obs_jj_temp_b(j,:)=squeeze(obs_temp_b(sp_jj(j),:));
    obs_jj_salt(j,:)=squeeze(obs_salt(sp_jj(j),:));
    obs_jj_salt_b(j,:)=squeeze(obs_salt_b(sp_jj(j),:));
    obs_jj_do(j,:)=squeeze(obs_do(sp_jj(j),:));
    obs_jj_do_b(j,:)=squeeze(obs_do_b(sp_jj(j),:));  
end


for i = 1:365
    %std
    obs_std_gy_no3(i)=nanstd(squeeze(obs_gy_no3(:,i)));
    obs_std_gy_no3_b(i)=nanstd(squeeze(obs_gy_no3_b(:,i)));
    obs_std_gy_nh4(i)=nanstd(squeeze(obs_gy_nh4(:,i)));
    obs_std_gy_nh4_b(i)=nanstd(squeeze(obs_gy_nh4_b(:,i)));
    obs_std_gy_chl(i)=nanstd(squeeze(obs_gy_chl(:,i)));
    obs_std_gy_chl_b(i)=nanstd(squeeze(obs_gy_chl_b(:,i)));
    obs_std_gy_temp(i)=nanstd(squeeze(obs_gy_temp(:,i)));
    obs_std_gy_temp_b(i)=nanstd(squeeze(obs_gy_temp_b(:,i)));
    obs_std_gy_salt(i)=nanstd(squeeze(obs_gy_salt(:,i)));
    obs_std_gy_salt_b(i)=nanstd(squeeze(obs_gy_salt_b(:,i)));
    obs_std_gy_do(i)=nanstd(squeeze(obs_gy_do(:,i)));
    obs_std_gy_do_b(i)=nanstd(squeeze(obs_gy_do_b(:,i)));

    obs_std_sgy_no3(i)=nanstd(squeeze(obs_sgy_no3(:,i)));
    obs_std_sgy_no3_b(i)=nanstd(squeeze(obs_sgy_no3_b(:,i)));
    obs_std_sgy_nh4(i)=nanstd(squeeze(obs_sgy_nh4(:,i)));
    obs_std_sgy_nh4_b(i)=nanstd(squeeze(obs_sgy_nh4_b(:,i)));
    obs_std_sgy_chl(i)=nanstd(squeeze(obs_sgy_chl(:,i)));
    obs_std_sgy_chl_b(i)=nanstd(squeeze(obs_sgy_chl_b(:,i)));
    obs_std_sgy_temp(i)=nanstd(squeeze(obs_sgy_temp(:,i)));
    obs_std_sgy_temp_b(i)=nanstd(squeeze(obs_sgy_temp_b(:,i)));
    obs_std_sgy_salt(i)=nanstd(squeeze(obs_sgy_salt(:,i)));
    obs_std_sgy_salt_b(i)=nanstd(squeeze(obs_sgy_salt_b(:,i)));
    obs_std_sgy_do(i)=nanstd(squeeze(obs_sgy_do(:,i)));
    obs_std_sgy_do_b(i)=nanstd(squeeze(obs_sgy_do_b(:,i)));

    obs_std_egy_no3(i)=nanstd(squeeze(obs_egy_no3(:,i)));
    obs_std_egy_no3_b(i)=nanstd(squeeze(obs_egy_no3_b(:,i)));
    obs_std_egy_nh4(i)=nanstd(squeeze(obs_egy_nh4(:,i)));
    obs_std_egy_nh4_b(i)=nanstd(squeeze(obs_egy_nh4_b(:,i)));
    obs_std_egy_chl(i)=nanstd(squeeze(obs_egy_chl(:,i)));
    obs_std_egy_chl_b(i)=nanstd(squeeze(obs_egy_chl_b(:,i)));
    obs_std_egy_temp(i)=nanstd(squeeze(obs_egy_temp(:,i)));
    obs_std_egy_temp_b(i)=nanstd(squeeze(obs_egy_temp_b(:,i)));
    obs_std_egy_salt(i)=nanstd(squeeze(obs_egy_salt(:,i)));
    obs_std_egy_salt_b(i)=nanstd(squeeze(obs_egy_salt_b(:,i)));
    obs_std_egy_do(i)=nanstd(squeeze(obs_egy_do(:,i)));
    obs_std_egy_do_b(i)=nanstd(squeeze(obs_egy_do_b(:,i)));

    obs_std_jj_no3(i)=nanstd(squeeze(obs_jj_no3(:,i)));
    obs_std_jj_no3_b(i)=nanstd(squeeze(obs_jj_no3_b(:,i)));
    obs_std_jj_nh4(i)=nanstd(squeeze(obs_jj_nh4(:,i)));
    obs_std_jj_nh4_b(i)=nanstd(squeeze(obs_jj_nh4_b(:,i)));
    obs_std_jj_chl(i)=nanstd(squeeze(obs_jj_chl(:,i)));
    obs_std_jj_chl_b(i)=nanstd(squeeze(obs_jj_chl_b(:,i)));
    obs_std_jj_temp(i)=nanstd(squeeze(obs_jj_temp(:,i)));
    obs_std_jj_temp_b(i)=nanstd(squeeze(obs_jj_temp_b(:,i)));
    obs_std_jj_salt(i)=nanstd(squeeze(obs_jj_salt(:,i)));
    obs_std_jj_salt_b(i)=nanstd(squeeze(obs_jj_salt_b(:,i)));
    obs_std_jj_do(i)=nanstd(squeeze(obs_jj_do(:,i)));
    obs_std_jj_do_b(i)=nanstd(squeeze(obs_jj_do_b(:,i)));
    
    %mean
    obm_gy_no3(i)=nanmean(squeeze(obs_gy_no3(:,i)));
    obm_gy_no3_b(i)=nanmean(squeeze(obs_gy_no3_b(:,i)));
    obm_gy_nh4(i)=nanmean(squeeze(obs_gy_nh4(:,i)));
    obm_gy_nh4_b(i)=nanmean(squeeze(obs_gy_nh4_b(:,i)));
    obm_gy_chl(i)=nanmean(squeeze(obs_gy_chl(:,i)));
    obm_gy_chl_b(i)=nanmean(squeeze(obs_gy_chl_b(:,i)));
    obm_gy_temp(i)=nanmean(squeeze(obs_gy_temp(:,i)));
    obm_gy_temp_b(i)=nanmean(squeeze(obs_gy_temp_b(:,i)));
    obm_gy_salt(i)=nanmean(squeeze(obs_gy_salt(:,i)));
    obm_gy_salt_b(i)=nanmean(squeeze(obs_gy_salt_b(:,i)));
    obm_gy_do(i)=nanmean(squeeze(obs_gy_do(:,i)));
    obm_gy_do_b(i)=nanmean(squeeze(obs_gy_do_b(:,i)));

    obm_sgy_no3(i)=nanmean(squeeze(obs_sgy_no3(:,i)));
    obm_sgy_no3_b(i)=nanmean(squeeze(obs_sgy_no3_b(:,i)));
    obm_sgy_nh4(i)=nanmean(squeeze(obs_sgy_nh4(:,i)));
    obm_sgy_nh4_b(i)=nanmean(squeeze(obs_sgy_nh4_b(:,i)));
    obm_sgy_chl(i)=nanmean(squeeze(obs_sgy_chl(:,i)));
    obm_sgy_chl_b(i)=nanmean(squeeze(obs_sgy_chl_b(:,i)));
    obm_sgy_temp(i)=nanmean(squeeze(obs_sgy_temp(:,i)));
    obm_sgy_temp_b(i)=nanmean(squeeze(obs_sgy_temp_b(:,i)));
    obm_sgy_salt(i)=nanmean(squeeze(obs_sgy_salt(:,i)));
    obm_sgy_salt_b(i)=nanmean(squeeze(obs_sgy_salt_b(:,i)));
    obm_sgy_do(i)=nanmean(squeeze(obs_sgy_do(:,i)));
    obm_sgy_do_b(i)=nanmean(squeeze(obs_sgy_do_b(:,i)));

    obm_egy_no3(i)=nanmean(squeeze(obs_egy_no3(:,i)));
    obm_egy_no3_b(i)=nanmean(squeeze(obs_egy_no3_b(:,i)));
    obm_egy_nh4(i)=nanmean(squeeze(obs_egy_nh4(:,i)));
    obm_egy_nh4_b(i)=nanmean(squeeze(obs_egy_nh4_b(:,i)));
    obm_egy_chl(i)=nanmean(squeeze(obs_egy_chl(:,i)));
    obm_egy_chl_b(i)=nanmean(squeeze(obs_egy_chl_b(:,i)));
    obm_egy_temp(i)=nanmean(squeeze(obs_egy_temp(:,i)));
    obm_egy_temp_b(i)=nanmean(squeeze(obs_egy_temp_b(:,i)));
    obm_egy_salt(i)=nanmean(squeeze(obs_egy_salt(:,i)));
    obm_egy_salt_b(i)=nanmean(squeeze(obs_egy_salt_b(:,i)));
    obm_egy_do(i)=nanmean(squeeze(obs_egy_do(:,i)));
    obm_egy_do_b(i)=nanmean(squeeze(obs_egy_do_b(:,i)));

    obm_jj_no3(i)=nanmean(squeeze(obs_jj_no3(:,i)));
    obm_jj_no3_b(i)=nanmean(squeeze(obs_jj_no3_b(:,i)));
    obm_jj_nh4(i)=nanmean(squeeze(obs_jj_nh4(:,i)));
    obm_jj_nh4_b(i)=nanmean(squeeze(obs_jj_nh4_b(:,i)));
    obm_jj_chl(i)=nanmean(squeeze(obs_jj_chl(:,i)));
    obm_jj_chl_b(i)=nanmean(squeeze(obs_jj_chl_b(:,i)));
    obm_jj_temp(i)=nanmean(squeeze(obs_jj_temp(:,i)));
    obm_jj_temp_b(i)=nanmean(squeeze(obs_jj_temp_b(:,i)));
    obm_jj_salt(i)=nanmean(squeeze(obs_jj_salt(:,i)));
    obm_jj_salt_b(i)=nanmean(squeeze(obs_jj_salt_b(:,i)));
    obm_jj_do(i)=nanmean(squeeze(obs_jj_do(:,i)));
    obm_jj_do_b(i)=nanmean(squeeze(obs_jj_do_b(:,i)));

end

%%plot

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obm_gy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3),obm_gy_no3./14 + obs_std_gy_no3./14,'m-','linew',2);
        plot(1:length(obm_gy_no3),obm_gy_no3./14 - obs_std_gy_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);
        plot(obm_gy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4),obm_gy_nh4./14 + obs_std_gy_nh4./14,'m-','linew',2);
        plot(1:length(obm_gy_nh4),obm_gy_nh4./14 - obs_std_gy_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);
        plot(obm_gy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
        plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);
        plot(obm_gy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp),obm_gy_temp + obs_std_gy_temp,'m-','linew',2);
        plot(1:length(obm_gy_temp),obm_gy_temp - obs_std_gy_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        plot(obm_gy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt),obm_gy_salt + obs_std_gy_salt,'m-','linew',2);
        plot(1:length(obm_gy_salt),obm_gy_salt - obs_std_gy_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        plot(obm_gy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 + obs_std_gy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 - obs_std_gy_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
        plot(obm_gy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 + obs_std_gy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 - obs_std_gy_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        plot(obm_gy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl_b),obm_gy_chl_b + obs_std_gy_chl_b,'m-','linew',2);
        plot(1:length(obm_gy_chl_b),obm_gy_chl_b - obs_std_gy_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);
        plot(obm_gy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp_b),obm_gy_temp_b + obs_std_gy_temp_b,'m-','linew',2);
        plot(1:length(obm_gy_temp_b),obm_gy_temp_b - obs_std_gy_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
        plot(obm_gy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt_b),obm_gy_salt_b + obs_std_gy_salt_b,'m-','linew',2);
        plot(1:length(obm_gy_salt_b),obm_gy_salt_b - obs_std_gy_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    
        
        
        
        

%plot daily
for i = 1:65
%   for   i = 4:4
   if ind_out(i) == 1
       fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obs_gy_no3(i,:)./14,'r-','linew',2); xlim([1 365]);
        title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        name_gf= name_tag{i};
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
        hold off;
        close;
   end
end

%plot daily
for i = 1:65
%   for   i = 4:4
   if ind_out(i) == 1
       fig = figure; hold on;
        plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(squeeze(mod_nh4(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
        plot(obs_nh4(i,:)./14,'r-','linew',2); xlim([1 365]);
        title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        name_gf= name_tag{i};
        print(fig,strcat('2001_daily_NH4_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
        hold off; close;
   end
end

%plot daily
for i = 1:65
%   for   i = 4:4
   if ind_out(i) == 1
       fig = figure; hold on;
        plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(squeeze(mod_do(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
        plot(obs_do(i,:).*0.7.*44.661,'r-','linew',2); xlim([1 365]);
        title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL DO']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('DO (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        name_gf= name_tag{i};
        print(fig,strcat('2001_daily_DO_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
        hold off;close;
   end
end

%plot daily
for i = 1:65
%   for   i = 4:4
   if ind_out(i) == 1
       fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(squeeze(mod_temp(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
        plot(obs_temp(i,:),'r-','linew',2); xlim([1 365]);
        title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('temp (^oC)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        name_gf= name_tag{i};
        print(fig,strcat('2001_daily_temp_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
        hold off;close;
   end
end

for i = 1:65
%   for   i = 4:4
   if ind_out(i) == 1
       fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(squeeze(mod_salt(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
        plot(obs_salt(i,:),'r-','linew',2); xlim([1 365]);
        title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        name_gf= name_tag{i};
        print(fig,strcat('2001_daily_salt_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
        hold off;close;
   end
end

%plot daily
for i = 1:65
%   for   i = 4:4
   if ind_out(i) == 1
       fig = figure; hold on;
        plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(squeeze(mod_chl(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
        plot(obs_chl(i,:),'r-','linew',2); xlim([1 365]);
        title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL chla']);
        xlabel('time(days on 2001)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        name_gf= name_tag{i};
        print(fig,strcat('2001_daily_chl_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
        hold off;
        close;
   end
end