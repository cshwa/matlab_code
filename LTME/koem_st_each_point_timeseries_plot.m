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

% make monthly mean
    
% obs_tdx = 49:60; %2001 %year
% obs_tdx = 61:72; %2002 %year
for j=1:65 % num. of st.

       obs_no3(j,:) = regime_no3(j,:);
       obs_nh4(j,:) = regime_nh4(j,:);
       obs_do(j,:) = regime_do(j,:);
       obs_chl(j,:) = regime_chl(j,:);
       obs_temp(j,:) = regime_temp(j,:);
       obs_salt(j,:) = regime_salt(j,:);
       % bot
       obs_no3_b(j,:) = regime_no3_b(j,:);
       obs_nh4_b(j,:) = regime_nh4_b(j,:);
       obs_do_b(j,:) = regime_do_b(j,:);
       obs_chl_b(j,:) = regime_chl_b(j,:);
       obs_temp_b(j,:) = regime_temp_b(j,:);
       obs_salt_b(j,:) = regime_salt_b(j,:);
end

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


for i = 1:length(ref_date)
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


return

figure; hold on;
plot(obs_gy_no3(end-2,:)./14,'ro');
plot(obs_gy_no3(end-1,:)./14,'bo');plot(obs_gy_no3(end,:)./14,'go')

plot(interp1(find(isnan(obs_gy_no3(end-2,:))==0),obs_gy_no3(end-2,find(isnan(obs_gy_no3(end-2,:))==0)),1:length(obs_gy_no3(end-2,:)))./14,'r')
hold on; plot(interp1(find(isnan(obs_gy_no3(end-1,:))==0),obs_gy_no3(end-1,find(isnan(obs_gy_no3(end-1,:))==0)),1:length(obs_gy_no3(end-1,:)))./14,'b')
hold on; plot(interp1(find(isnan(obs_gy_no3(end,:))==0),obs_gy_no3(end,find(isnan(obs_gy_no3(end,:))==0)),1:length(obs_gy_no3(end,:)))./14,'g')
xlim([1 length(obs_gy_no3(end,:))])


figure;
for i = 1:size(obs_gy_no3,1)
plot(interp1(find(isnan(obs_gy_no3(i,:))==0),obs_gy_no3(i,find(isnan(obs_gy_no3(i,:))==0)),1:length(obs_gy_no3(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_no3(end,:))])



figure;
for i = 1:size(obs_gy_nh4,1)
plot(interp1(find(isnan(obs_gy_nh4(i,:))==0),obs_gy_nh4(i,find(isnan(obs_gy_nh4(i,:))==0)),1:length(obs_gy_nh4(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_nh4(end,:))])


figure;
for i = 1:size(obs_gy_chl,1)
plot(interp1(find(isnan(obs_gy_chl(i,:))==0),obs_gy_chl(i,find(isnan(obs_gy_chl(i,:))==0)),1:length(obs_gy_chl(i,:))),'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_chl(end,:))])

% 3 point
% extract 2, 5, 8, 11 month
obs_gy_no3_f=NaN(size(obs_gy_no3,1),length(obs_gy_no3(1,2:12:end))*4);
obs_gy_no3_f(:,1:4:end)=obs_gy_no3(:,2:12:end);
obs_gy_no3_f(:,2:4:end)=obs_gy_no3(:,5:12:end);
obs_gy_no3_f(:,3:4:end)=obs_gy_no3(:,8:12:end);
obs_gy_no3_f(:,4:4:end)=obs_gy_no3(:,11:12:end);

obs_gy_nh4_f=NaN(size(obs_gy_nh4,1),length(obs_gy_nh4(1,2:12:end))*4);
obs_gy_nh4_f(:,1:4:end)=obs_gy_nh4(:,2:12:end);
obs_gy_nh4_f(:,2:4:end)=obs_gy_nh4(:,5:12:end);
obs_gy_nh4_f(:,3:4:end)=obs_gy_nh4(:,8:12:end);
obs_gy_nh4_f(:,4:4:end)=obs_gy_nh4(:,11:12:end);

obs_gy_chl_f=NaN(size(obs_gy_chl,1),length(obs_gy_chl(1,2:12:end))*4);
obs_gy_chl_f(:,1:4:end)=obs_gy_chl(:,2:12:end);
obs_gy_chl_f(:,2:4:end)=obs_gy_chl(:,5:12:end);
obs_gy_chl_f(:,3:4:end)=obs_gy_chl(:,8:12:end);
obs_gy_chl_f(:,4:4:end)=obs_gy_chl(:,11:12:end);


figure;
for i = 7:size(obs_gy_no3_f,1)
plot(interp1(find(isnan(obs_gy_no3_f(i,:))==0),obs_gy_no3_f(i,find(isnan(obs_gy_no3_f(i,:))==0)),1:length(obs_gy_no3_f(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_no3_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS no3']);
        xticks(1:4:length(obs_gy_no3_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('no3 (umol/m^3)','fontsize',13)
        grid on
        legend(name_tag{33:35})
        ylim([0 20])
        set(gca,'fontsize',13)
        xtickangle(60)

figure;
for i = 7:size(obs_gy_nh4_f,1)
plot(interp1(find(isnan(obs_gy_nh4_f(i,:))==0),obs_gy_nh4_f(i,find(isnan(obs_gy_nh4_f(i,:))==0)),1:length(obs_gy_nh4_f(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_nh4_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS nh4']);
        xticks(1:4:length(obs_gy_nh4_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 10])
        legend(name_tag{33:35})
        set(gca,'fontsize',13)
        xtickangle(60)
        
figure;
for i = 7:size(obs_gy_chl_f,1)
plot(interp1(find(isnan(obs_gy_chl_f(i,:))==0),obs_gy_chl_f(i,find(isnan(obs_gy_chl_f(i,:))==0)),1:length(obs_gy_chl_f(i,:))),'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_chl_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS chl']);
        xticks(1:4:length(obs_gy_chl_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('chl (umol/m^3)','fontsize',13)
        grid on
        ylim([0 20])
        legend(name_tag{33:35})
        set(gca,'fontsize',13)
        xtickangle(60)

%% all sp_gy point
        clearvars na_gy
        for i = 1:length(sp_gy)
        na_gy{i} = string(name_tag{sp_gy(i)});
        end
        
figure;
for i = 1:size(obs_gy_no3_f,1)
plot(interp1(find(isnan(obs_gy_no3_f(i,:))==0),obs_gy_no3_f(i,find(isnan(obs_gy_no3_f(i,:))==0)),1:length(obs_gy_no3_f(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_no3_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS no3']);
        xticks(1:4:length(obs_gy_no3_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('no3 (umol/m^3)','fontsize',13)
        grid on
        legend(na_gy)
        ylim([0 35])
        set(gca,'fontsize',13)
        xtickangle(60)

figure;
for i = 1:size(obs_gy_nh4_f,1)
plot(interp1(find(isnan(obs_gy_nh4_f(i,:))==0),obs_gy_nh4_f(i,find(isnan(obs_gy_nh4_f(i,:))==0)),1:length(obs_gy_nh4_f(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_nh4_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS nh4']);
        xticks(1:4:length(obs_gy_nh4_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 12])
        legend(na_gy)
        set(gca,'fontsize',13)
        xtickangle(60)
        
figure;
for i = 1:size(obs_gy_chl_f,1)
plot(interp1(find(isnan(obs_gy_chl_f(i,:))==0),obs_gy_chl_f(i,find(isnan(obs_gy_chl_f(i,:))==0)),1:length(obs_gy_chl_f(i,:))),'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_chl_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS chl']);
        xticks(1:4:length(obs_gy_chl_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('chl (umol/m^3)','fontsize',13)
        grid on
        ylim([0 20])
        legend(na_gy)
        set(gca,'fontsize',13)
        xtickangle(60)
%% 4 point
sp_pic=[1,3,4,5];
name_tag{sp_gy(sp_pic)};

clearvars na_gy
for i = 1:4
na_gy{i} = string(name_tag{sp_gy(sp_pic(i))});
end

figure;
for i = 1:4
plot(interp1(find(isnan(obs_gy_no3_f(sp_pic(i),:))==0),obs_gy_no3_f(sp_pic(i),find(isnan(obs_gy_no3_f(sp_pic(i),:))==0)),1:length(obs_gy_no3_f(sp_pic(i),:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_no3_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS no3']);
        xticks(1:4:length(obs_gy_no3_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('no3 (umol/m^3)','fontsize',13)
        grid on
        legend(na_gy)
        ylim([0 35])
        set(gca,'fontsize',13)
        xtickangle(60)

figure;
for i = 1:4
plot(interp1(find(isnan(obs_gy_nh4_f(sp_pic(i),:))==0),obs_gy_nh4_f(sp_pic(i),find(isnan(obs_gy_nh4_f(sp_pic(i),:))==0)),1:length(obs_gy_nh4_f(sp_pic(i),:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_nh4_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS nh4']);
        xticks(1:4:length(obs_gy_nh4_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 10])
        legend(na_gy)
        set(gca,'fontsize',13)
        xtickangle(60)
        
figure;
for i = 1:4
plot(interp1(find(isnan(obs_gy_chl_f(sp_pic(i),:))==0),obs_gy_chl_f(sp_pic(i),find(isnan(obs_gy_chl_f(sp_pic(i),:))==0)),1:length(obs_gy_chl_f(sp_pic(i),:))),'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_chl_f(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS chl']);
        xticks(1:4:length(obs_gy_chl_f))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('chl (umol/m^3)','fontsize',13)
        grid on
        ylim([0 20])
        legend(na_gy)
        set(gca,'fontsize',13)
        xtickangle(60)


% 3point old
figure;
for i = 7:size(obs_gy_no3,1)
plot(interp1(find(isnan(obs_gy_no3(i,:))==0),obs_gy_no3(i,find(isnan(obs_gy_no3(i,:))==0)),1:length(obs_gy_no3(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_no3(end,:))])
clearvars ref_plt_t
j=0
for i = 1:12:length(ref_date)
    j=j+1;
    ref_plt_t{j}=ref_date{i}(1:4);
end    
        title(['KOEM OBS no3']);
        xticks(1:12:length(ref_date))
        xticklabels(ref_plt_t)
        xlabel('time','fontsize',13)
        ylabel('no3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 20])
        set(gca,'fontsize',13)
        xtickangle(90)

figure;
for i = 7:size(obs_gy_nh4,1)
plot(interp1(find(isnan(obs_gy_nh4(i,:))==0),obs_gy_nh4(i,find(isnan(obs_gy_nh4(i,:))==0)),1:length(obs_gy_nh4(i,:)))./14,'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_nh4(end,:))])


figure;
for i = 7:size(obs_gy_chl,1)
plot(interp1(find(isnan(obs_gy_chl(i,:))==0),obs_gy_chl(i,find(isnan(obs_gy_chl(i,:))==0)),1:length(obs_gy_chl(i,:))),'linew',2)
hold on;
end
alpha(0.3)
xlim([1 length(obs_gy_chl(end,:))])