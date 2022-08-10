close all; clear; clc; 
k = [2012,2013];

% % cd G:\장기생태\Dynamic\result\2011\Out\monthly
% for i = 1:7
%      for j =  1:12
%         h=k(i); 
%         temp(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\monthly\monthly_',num2str(h,'%04d'),'_',num2str(j,'%02d'),'.nc'],'temp');
%         salt(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\monthly\monthly_',num2str(h,'%04d'),'_',num2str(j,'%02d'),'.nc'],'salt');
%         zeta(:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\monthly\monthly_',num2str(h,'%04d'),'_',num2str(j,'%02d'),'.nc'],'zeta');
%         u(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\monthly\monthly_',num2str(h,'%04d'),'_',num2str(j,'%02d'),'.nc'],'u');
%         v(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\monthly\monthly_',num2str(h,'%04d'),'_',num2str(j,'%02d'),'.nc'],'v');
%      end
% end
% 

temp = NaN(252,176,20,length(k),366);
salt = NaN(252,176,20,length(k),366);
zeta = NaN(252,176,length(k),366);

leapyr = leapyear(k);

for i = 1:length(k)
     if leapyr(i) == 1
         for j = 1:366
            h=k(i); 
            temp(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'temp');
            salt(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'salt');
            zeta(:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'zeta'); 
         end
     else
         for j =  1:365
            h=k(i); 
            temp(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'temp');
            salt(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'salt');
            zeta(:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'zeta');
         end
     end
end
for i = 1:366
    temp_366(:,:,:,i)=nanmean(squeeze(temp(:,:,:,:,i)),4);
    salt_366(:,:,:,i)=nanmean(squeeze(salt(:,:,:,:,i)),4);
    zeta_366(:,:,i)=nanmean(squeeze(zeta(:,:,:,i)),3);
end

save('2012_2013_TS_mean_GY.mat','temp_366','salt_366','-v7.3');
save('2012_2013_TS_GY.mat','-v7.3');

close all; clear; clc; 
k = [2012,2013];

leapyr = leapyear(k);

u = NaN(251,176,20,length(k),366);
v = NaN(252,175,20,length(k),366);
for i = 1:length(k)
         if leapyr(i) == 1
            for j = 1:366
             h=k(i);
             u(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'u');
             v(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'v');
            end
         else
            for j =  1:365
            h=k(i);
            u(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'u');
            v(:,:,:,i,j) = ncread(['G:\장기생태\NPZD\',num2str(h,'%04d'),'\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'v');
            end
         end
end
 
 for i = 1:366
    u_366(:,:,:,i)=nanmean(squeeze(u(:,:,:,:,i)),4);
    v_366(:,:,:,i)=nanmean(squeeze(v(:,:,:,:,i)),4);
 end

 
save('2012_2013_uv_mean_GY.mat','u_366','v_366','-v7.3');
save('2012_2013_uv_GY.mat','-v7.3');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc
% make 1980~present
load 2012_2013_uv_mean_GY.mat
k=0
for i = 1980:1980
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

clearvars em em_seanson k_f kk
k_f=[1:12];
for i = 1:length(k_f)
    em{i} = sum(eom_d(1:k_f(i)-1))+1:sum(eom_d(1:k_f(i)));
end

for i = 1:12
    u_mth_pre(:,:,:,i)=nanmean(u_366(:,:,:,em{i}),4);
    v_mth_pre(:,:,:,i)=nanmean(v_366(:,:,:,em{i}),4);
end

clearvars u v
k = [2014,2015,2016,2017];
for i = 1 : length(k)
for j =  1:12
        u(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'u');
        v(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'v');
end
end

u_mth_pre2 = squeeze(nanmean(u,5));
v_mth_pre2 = squeeze(nanmean(v,5));

u_mth= (u_mth_pre + u_mth_pre2)./2;
v_mth= (v_mth_pre + v_mth_pre2)./2;

save('2012_2017_uv_mean_mth_GY.mat','*_mth','-v7.3');
save('2012_2017_uv_GY_mth.mat','-v7.3');


close all; clear; clc

load 2012_2013_TS_GY.mat

% make 1980~present
k=0
for i = 1980:1980
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

clearvars em em_seanson k_f kk
k_f=[1:12];
for i = 1:length(k_f)
    em{i} = sum(eom_d(1:k_f(i)-1))+1:sum(eom_d(1:k_f(i)));
end

for i = 1:12
    temp_mth_pre(:,:,:,i)=nanmean(temp_366(:,:,:,em{i}),4);
    salt_mth_pre(:,:,:,i)=nanmean(salt_366(:,:,:,em{i}),4);
    zeta_mth_pre(:,:,i)=nanmean(zeta_366(:,:,em{i}),3);
end

clearvars temp salt zeta
k = [2014,2015,2016,2017];
for i = 1 : length(k)
for j =  1:12
        temp(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'temp');
        salt(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'salt');
        zeta(:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'zeta');
end
end

temp_mth_pre2 = squeeze(nanmean(temp,5));
salt_mth_pre2 = squeeze(nanmean(salt,5));
zeta_mth_pre2 = squeeze(nanmean(zeta,4));
 
temp_mth= (temp_mth_pre + temp_mth_pre2)./2;
salt_mth= (salt_mth_pre + salt_mth_pre2)./2;
zeta_mth= (zeta_mth_pre + zeta_mth_pre2)./2;

save('2012_2017_TS_mean_mth_GY.mat','*_mth','-v7.3');
save('2012_2017_TS_GY_mth.mat','-v7.3');


%% bio

close all; clear; clc

k = [2015,2016,2017];
for i = 1 : length(k)
for j =  1:12
        no3(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'NO3');
        nh4(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'NH4');
end
end

no3_mth = squeeze(nanmean(no3,5));
nh4_mth = squeeze(nanmean(nh4,5));

save('2015_2017_nitro_mean_mth_GY.mat','*_mth','-v7.3');
save('2015_2017_nitro_GY_mth.mat','-v7.3');


close all; clear; clc

k = [2015,2016,2017];
for i = 1 : length(k)
for j =  1:12
        chla(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'chlorophyll');
        phy(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'phytoplankton');
        zoo(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'zooplankton');
        do(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'oxygen');
end
end

chla_mth = squeeze(nanmean(chla,5));
phy_mth = squeeze(nanmean(phy,5));
zoo_mth = squeeze(nanmean(zoo,5));
do_mth = squeeze(nanmean(do,5));

save('2015_2017_bio_mean_mth_GY.mat','*_mth','-v7.3');
save('2015_2017_bio_GY_mth.mat','-v7.3');

close all; clear; clc

k = [2015,2016,2017];
for i = 1 : length(k)
for j =  1:12
        ld(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'LdetritusN');
        sd(:,:,:,j,i) = ncread(['G:\장기생태\NPZD\' num2str(k(i),'%04d') '\monthly\monthly_' num2str(k(i),'%04d') '_',num2str(j,'%02d'),'.nc'],'SdetritusN');
end
end

ld_mth = squeeze(nanmean(ld,5));
sd_mth = squeeze(nanmean(sd,5));

save('2015_2017_bio2_mean_mth_GY.mat','*_mth','-v7.3');
save('2015_2017_bio2_GY_mth.mat','-v7.3');





