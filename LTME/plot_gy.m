%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
k = [1980:1986];
for i = 1:7
     for j =  1:365
        h=k(i); 
        temp(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/',num2str(h,'%04d'),'_regrid_v3/ocean_avg_',num2str(j,'%04d'),'.nc'],'temp');
        salt(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/',num2str(h,'%04d'),'_regrid_v3/ocean_avg_',num2str(j,'%04d'),'.nc'],'salt');
        zeta(:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/',num2str(h,'%04d'),'_regrid_v3/ocean_avg_',num2str(j,'%04d'),'.nc'],'zeta');
     end
end


for i = 1:7
     for j =  1:365
        h=k(i); 
        u(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/',num2str(h,'%04d'),'_regrid_v3/ocean_avg_',num2str(j,'%04d'),'.nc'],'u');
        v(:,:,:,i,j) = ncread(['/data1/cshwa/ext_hdd2/longterm/coco_storage/',num2str(h,'%04d'),'_regrid_v3/ocean_avg_',num2str(j,'%04d'),'.nc'],'v');
     end
end

for i = 1:365
    temp_365(:,:,:,i)=mean(squeeze(temp(:,:,:,:,i)),4);
    salt_365(:,:,:,i)=mean(squeeze(salt(:,:,:,:,i)),4);
    u_365(:,:,:,i)=mean(squeeze(u(:,:,:,:,i)),4);
    v_365(:,:,:,i)=mean(squeeze(v(:,:,:,:,i)),4);
    zeta_365(:,:,:,i)=mean(squeeze(zeta(:,:,:,i)),3);
end

save('1980_1986_mean_GY.mat','-v7.3');

save('1980_1986_GY.mat','temp_365','salt_365','u_365','v_365','-v7.3');


% make 1980~present
k=0
for i = 1981:1981
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

clearvars *_mth
for i = 1:12
    temp_mth(:,:,:,i)=nanmean(temp_365(:,:,:,em{i}),4);
    salt_mth(:,:,:,i)=nanmean(salt_365(:,:,:,em{i}),4);
    u_mth(:,:,:,i)=nanmean(u_365(:,:,:,em{i}),4);
    v_mth(:,:,:,i)=nanmean(v_365(:,:,:,em{i}),4);
    zeta_mth(:,:,i)=nanmean(squeeze(zeta_365(:,:,1,em{i})),3);
end

save('1980_1986_mean_mth_all_GY.mat','-v7.3');
save('1980_1986_mean_mth_GY.mat','temp_mth','salt_mth','u_mth','v_mth','-v7.3');








