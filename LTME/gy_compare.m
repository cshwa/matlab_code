close all; clear; clc; 
k = [2011,2012,2013,2015,2016];

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

temp = NaN(252,176,20,length(k),365);
salt = NaN(252,176,20,length(k),365);
zeta = NaN(252,176,length(k),365);


for i = 1:length(k)
     for j =  1:365
        h=k(i); 
        temp(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'temp');
        salt(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'salt');
        zeta(:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'zeta');
     end
end
for i = 1:365
    temp_365(:,:,:,i)=nanmean(squeeze(temp(:,:,:,:,i)),4);
    salt_365(:,:,:,i)=nanmean(squeeze(salt(:,:,:,:,i)),4);
    zeta_365(:,:,i)=nanmean(squeeze(zeta(:,:,:,i)),3);
end

save('2011_2016_TS_mean_GY.mat','temp_365','salt_365','-v7.3');
save('2011_2016_TS_GY.mat','-v7.3');

close all; clear; clc; 
k = [2011,2012,2013,2015,2016];

u = NaN(251,176,20,length(k),365);
v = NaN(252,175,20,length(k),365);
for i = 1:length(k)
     for j =  1:365
       h=k(i);
 u(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'u');
 v(:,:,:,i,j) = ncread(['G:\장기생태\Dynamic\result\',num2str(h,'%04d'),'\Out\model_result\ocean_avg_',num2str(j,'%04d'),'.nc'],'v');
     end
end
 
 for i = 1:365
    u_365(:,:,:,i)=nanmean(squeeze(u(:,:,:,:,i)),4);
    v_365(:,:,:,i)=nanmean(squeeze(v(:,:,:,:,i)),4);
 end

 
save('2011_2016_uv_mean_GY.mat','u_365','v_365','-v7.3');
save('2011_2016_uv_GY.mat','-v7.3');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc
% make 1980~present
load 2011_2016_uv_mean_GY.mat
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

for i = 1:12
    u_mth_pre(:,:,:,i)=nanmean(u_365(:,:,:,em{i}),4);
    v_mth_pre(:,:,:,i)=nanmean(v_365(:,:,:,em{i}),4);
end

for j =  1:12
        u(:,:,:,j) = ncread(['G:\장기생태\Dynamic\result\2014\monthly\monthly_2014_',num2str(j,'%02d'),'.nc'],'u');
        v(:,:,:,j) = ncread(['G:\장기생태\Dynamic\result\2014\monthly\monthly_2014_',num2str(j,'%02d'),'.nc'],'v');
end


u_mth= (u_mth_pre + u)./2;
v_mth= (v_mth_pre + v)./2;

save('2011_2016_uv_mean_mth_GY.mat','*_mth','-v7.3');
save('2011_2016_uv_GY_mth.mat','-v7.3');


close all; clear; clc

load 2011_2016_TS_GY.mat

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

for i = 1:12
    temp_mth_pre(:,:,:,i)=nanmean(temp_365(:,:,:,em{i}),4);
    salt_mth_pre(:,:,:,i)=nanmean(salt_365(:,:,:,em{i}),4);
    zeta_mth_pre(:,:,i)=nanmean(zeta_365(:,:,em{i}),3);
end

clearvars temp salt zeta
for j =  1:12
        temp(:,:,:,j) = ncread(['G:\장기생태\Dynamic\result\2014\monthly\monthly_2014_',num2str(j,'%02d'),'.nc'],'temp');
        salt(:,:,:,j) = ncread(['G:\장기생태\Dynamic\result\2014\monthly\monthly_2014_',num2str(j,'%02d'),'.nc'],'salt');
        zeta(:,:,j) = ncread(['G:\장기생태\Dynamic\result\2014\monthly\monthly_2014_',num2str(j,'%02d'),'.nc'],'zeta');
end
 
temp_mth= (temp_mth_pre + temp)./2;
salt_mth= (salt_mth_pre + salt)./2;
zeta_mth= (zeta_mth_pre + zeta)./2;

save('2011_2016_TS_mean_mth_GY.mat','*_mth','-v7.3');
save('2011_2016_TS_GY_mth.mat','-v7.3');




