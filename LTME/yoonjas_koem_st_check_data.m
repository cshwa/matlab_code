close all; clc; clear;

% name_tag{sp_gy}
% ê´‘ì–‘?•­, ê´‘ì–‘4, ê´‘ì–‘3, ê´‘ì–‘2, ê´‘ì–‘1, ê´‘ì–‘5, ?—¬?ˆ˜2, ?—¬?ˆ˜3, ?—¬?ˆ˜1
% [res I]=sort([4,3,2,1,5]);
P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

yoonja=load('yoonjakangs_koem_data_monthly.mat');

[res I]=sort([1,5,4,3,2,6,8,9,7]);
% sorting confirm
nn=['±¤¾çÇ×'; '±¤¾ç4'; '±¤¾ç3'; '±¤¾ç2'; '±¤¾ç1'; '±¤¾ç5'; '¿©¼ö2'; '¿©¼ö3'; '¿©¼ö1'];
nn(I,:)

% % nn(I,:)
loc(1,1:2) = [34.8625;	127.7403];
loc(2,1:2) = [34.85139;  127.6797];
loc(3,1:2) = [34.88389;  127.6506];
loc(4,1:2) = [34.90222;  127.6822];
loc(5,1:2) = [34.92083;  127.8233];
loc(6,1:2) = [34.83194;  127.8011];
loc(7,1:2) = [34.73556;  127.7661];
loc(8,1:2) = [34.76278;  127.7614];
loc(9,1:2) = [34.76444;  127.8053];

grd_file='grid_gy_v11_s.nc';
x=ncread(grd_file, 'lon_rho');
y=ncread(grd_file, 'lat_rho');
mask=ncread(grd_file, 'mask_rho');
h=ncread(grd_file, 'h');

%plot only in the domain
figure; hold on; %pcolor(lon,lat,(mask./mask)); shading flat;
set(gca,'fontsize',15,'fontweight','bold'); 
grid on; 
for i=1:length(nn)
        text(loc(i,2)+0.005,loc(i,1), char(nn(I(i),:)),'color','r'); %plot only 97
        plot(loc(i,2),loc(i,1),'.','color','r'); %plot only 97
end
ylim([34.723 35])
xlim([min(min(x)) 127.9])
% xlim([min(min(lon)) max(max(lon))])
% ylim([min(min(lat)) max(max(lat))])
plot_google_map('MapType','terrain','Scale',2,'Resize',2)   % overlay google map
ylim([34.723 35])
xlim([min(min(x)) 127.9])

fig=figure; hold on;
for i = 1:4
temp_data = yoonja.nh4_sur(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 inf]);
ylabel('NH4-N (umol/L)')
title('NH4 surface')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('NH4_s_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = yoonja.nh4_sur(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 inf]);
ylabel('NH4-N (umol/L)')
title('NH4 surface')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('NH4_s_6to9.png'),'-dpng')

%% make year
for i = 1:9
    for j=1:length(1997:2019)
       nh4_s_yr(i,j) = nanmean(yoonja.nh4_sur(i,1+((j-1)*12):12*j));
       nh4_b_yr(i,j) = nanmean(yoonja.nh4_bot(i,1+((j-1)*12):12*j));
       no3_s_yr(i,j) = nanmean(yoonja.no3_sur(i,1+((j-1)*12):12*j));
       no3_b_yr(i,j) = nanmean(yoonja.no3_bot(i,1+((j-1)*12):12*j));
       chl_s_yr(i,j) = nanmean(yoonja.chl_sur(i,1+((j-1)*12):12*j));
       chl_b_yr(i,j) = nanmean(yoonja.chl_bot(i,1+((j-1)*12):12*j));
       po4_s_yr(i,j) = nanmean(yoonja.po4_sur(i,1+((j-1)*12):12*j));
       po4_b_yr(i,j) = nanmean(yoonja.po4_bot(i,1+((j-1)*12):12*j));
       temp_s_yr(i,j) = nanmean(yoonja.temp_sur(i,1+((j-1)*12):12*j));
       temp_b_yr(i,j) = nanmean(yoonja.temp_bot(i,1+((j-1)*12):12*j));
       salt_s_yr(i,j) = nanmean(yoonja.salt_sur(i,1+((j-1)*12):12*j));
       salt_b_yr(i,j) = nanmean(yoonja.salt_bot(i,1+((j-1)*12):12*j));
    end
end

for i = 1:9
    for j=1:length(1997:2019)
       nh4_s_mn(i,j) = nanmean(yoonja.nh4_sur(i,1+((j-1)*12):12*j));
       nh4_b_yr(i,j) = nanmean(yoonja.nh4_bot(i,1+((j-1)*12):12*j));
       no3_s_yr(i,j) = nanmean(yoonja.no3_sur(i,1+((j-1)*12):12*j));
       no3_b_yr(i,j) = nanmean(yoonja.no3_bot(i,1+((j-1)*12):12*j));
       chl_s_yr(i,j) = nanmean(yoonja.chl_sur(i,1+((j-1)*12):12*j));
       chl_b_yr(i,j) = nanmean(yoonja.chl_bot(i,1+((j-1)*12):12*j));
       po4_s_yr(i,j) = nanmean(yoonja.po4_sur(i,1+((j-1)*12):12*j));
       po4_b_yr(i,j) = nanmean(yoonja.po4_bot(i,1+((j-1)*12):12*j));
       temp_s_yr(i,j) = nanmean(yoonja.temp_sur(i,1+((j-1)*12):12*j));
       temp_b_yr(i,j) = nanmean(yoonja.temp_bot(i,1+((j-1)*12):12*j));
       salt_s_yr(i,j) = nanmean(yoonja.salt_sur(i,1+((j-1)*12):12*j));
       salt_b_yr(i,j) = nanmean(yoonja.salt_bot(i,1+((j-1)*12):12*j));
    end
end

%% surface

%% temp
fig=figure; hold on;
for i = 1:4
temp_data = temp_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('temp_s_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = temp_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('temp_s_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = temp_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('temp_s_6to9.png'),'-dpng')

%% salt
fig=figure; hold on;
for i = 1:4
salt_data = salt_s_yr(i,:);
no_nan_data= salt_data(~isnan(salt_data));
no_nan_ind=find(isnan(salt_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('salt_s_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
salt_data = salt_s_yr(i,:);
no_nan_data= salt_data(~isnan(salt_data));
no_nan_ind=find(isnan(salt_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('salt_s_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
salt_data = salt_s_yr(i,:);
no_nan_data= salt_data(~isnan(salt_data));
no_nan_ind=find(isnan(salt_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('salt_s_6to9.png'),'-dpng')

%% nh4
fig=figure; hold on;
for i = 1:4
temp_data = nh4_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 surface')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('NH4_s_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = nh4_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 surface')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('NH4_s_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = nh4_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 surface')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('NH4_s_6to9.png'),'-dpng')

%% no3
fig=figure; hold on;
for i = 1:4
temp_data = no3_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('no3-N (umol/L)')
title('no3 surface')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('no3_s_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = no3_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('no3-N (umol/L)')
title('no3 surface')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('no3_s_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = no3_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('no3-N (umol/L)')
title('no3 surface')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('no3_s_6to9.png'),'-dpng')


%% po4
fig=figure; hold on;
for i = 1:4
temp_data = po4_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('po4-N (umol/L)')
title('po4 surface')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('po4_s_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = po4_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('po4-N (umol/L)')
title('po4 surface')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('po4_s_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = po4_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('po4-N (umol/L)')
title('po4 surface')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('po4_s_6to9.png'),'-dpng')


%% chl
fig=figure; hold on;
for i = 1:4
temp_data = chl_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('chl_s_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = chl_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('chl_s_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = chl_s_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('chl_s_6to9.png'),'-dpng')


%% bottom
%% temp
fig=figure; hold on;
for i = 1:4
temp_data = temp_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('temp_b_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = temp_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('temp_b_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = temp_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('temp_b_6to9.png'),'-dpng')

%% salt
fig=figure; hold on;
for i = 1:4
salt_data = salt_b_yr(i,:);
no_nan_data= salt_data(~isnan(salt_data));
no_nan_ind=find(isnan(salt_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('salt_b_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
salt_data = salt_b_yr(i,:);
no_nan_data= salt_data(~isnan(salt_data));
no_nan_ind=find(isnan(salt_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('salt_b_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
salt_data = salt_b_yr(i,:);
no_nan_data= salt_data(~isnan(salt_data));
no_nan_ind=find(isnan(salt_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('salt_b_6to9.png'),'-dpng')

%% nh4
fig=figure; hold on;
for i = 1:4
temp_data = nh4_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 bottom')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('NH4_b_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = nh4_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 bottom')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('NH4_b_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = nh4_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 bottom')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('NH4_b_6to9.png'),'-dpng')

%% no3
fig=figure; hold on;
for i = 1:4
temp_data = no3_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('NO3-N (umol/L)')
title('no3 bottom')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('no3_b_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = no3_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('NO3-N (umol/L)')
title('no3 bottom')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('no3_b_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = no3_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('NO3-N (umol/L)')
title('no3 bottom')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('no3_b_6to9.png'),'-dpng')


%% po4
fig=figure; hold on;
for i = 1:4
temp_data = po4_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('PO4-P (umol/L)')
title('po4 bottom')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('po4_b_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = po4_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('PO4-P (umol/L)')
title('po4 bottom')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('po4_b_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = po4_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('PO4-P (umol/L)')
title('po4 bottom')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('po4_b_6to9.png'),'-dpng')


%% chl
fig=figure; hold on;
for i = 1:4
temp_data = chl_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
legend(nn(I(1:4),:))
print(fig,strcat('chl_b_1to4.png'),'-dpng')

fig=figure; hold on;
for i = 5
temp_data = chl_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
legend(nn(I(5),:))
print(fig,strcat('chl_b_5.png'),'-dpng')

fig=figure; hold on;
for i = 6:9
temp_data = chl_b_yr(i,:);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
end
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
legend(nn(I(6:9),:))
print(fig,strcat('chl_b_6to9.png'),'-dpng')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% space mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% surface

%% temp
fig=figure; hold on;
i = 1:4
temp_data = nanmean(temp_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(temp_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('temp_s_sp_mean.png'),'-dpng')


fig=figure; hold on;
i = 2:4
temp_data = nanmean(temp_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(temp_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('temp_s_sp_mean(minus1).png'),'-dpng')

%% salt
fig=figure; hold on;
i = 1:4
temp_data = nanmean(salt_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(salt_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('salt_s_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(salt_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(salt_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('salt_s_sp_mean(minus1).png'),'-dpng')


%% nh4
fig=figure; hold on;
i = 1:4
temp_data = nanmean(nh4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(nh4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('nh4_s_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(nh4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(nh4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('nh4_s_sp_mean(minus1).png'),'-dpng')


%% no3
fig=figure; hold on;
i = 1:4
temp_data = nanmean(no3_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('NO3-N (umol/L)')
title('NO3 surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(no3_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('no3_s_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(no3_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('NO3-N (umol/L)')
title('NO3 surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(no3_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('no3_s_sp_mean(minus1).png'),'-dpng')


%% po4
fig=figure; hold on;
i = 1:4
temp_data = nanmean(po4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('po4-P (umol/L)')
title('po4 surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(po4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('po4_s_sp_mean.png'),'-dpng')


fig=figure; hold on;
i = 2:4
temp_data = nanmean(po4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('po4-P (umol/L)')
title('po4 surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(po4_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('po4_s_sp_mean(minus1).png'),'-dpng')


%% chl
fig=figure; hold on;
i = 1:4
temp_data = nanmean(chl_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(chl_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('chl_s_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(chl_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl surface')
xtickangle(45);
i = 6:9
temp_data = nanmean(chl_s_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('chl_s_sp_mean(minus1).png'),'-dpng')


%% bottom
%% temp
fig=figure; hold on;
i = 1:4
temp_data = nanmean(temp_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(temp_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('temp_b_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(temp_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([14 20]);
ylabel('temp (^oC)')
title('temp bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(temp_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('temp_b_sp_mean(minus1).png'),'-dpng')


%% salt
fig=figure; hold on;
i = 1:4
temp_data = nanmean(salt_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(salt_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('salt_b_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(salt_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([26 34]);
ylabel('salt (psu)')
title('salt bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(salt_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('salt_b_sp_mean(minus1).png'),'-dpng')

%% nh4
fig=figure; hold on;
i = 1:4
temp_data = nanmean(nh4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(nh4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('nh4_b_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(nh4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('NH4-N (umol/L)')
title('NH4 bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(nh4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('nh4_b_sp_mean(minus1).png'),'-dpng')


%% no3
fig=figure; hold on;
i = 1:4
temp_data = nanmean(no3_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('NO3-N (umol/L)')
title('NO3 bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(no3_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('no3_b_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(no3_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 20]);
ylabel('NO3-N (umol/L)')
title('NO3 bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(no3_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./N_MW ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('no3_b_sp_mean(minus1).png'),'-dpng')


%% po4
fig=figure; hold on;
i = 1:4
temp_data = nanmean(po4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('po4-P (umol/L)')
title('po4 bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(po4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('po4_b_sp_mean.png'),'-dpng')


fig=figure; hold on;
i = 2:4
temp_data = nanmean(po4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 8]);
ylabel('po4-P (umol/L)')
title('po4 bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(po4_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data./P_MW ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('po4_b_sp_mean(minus1).png'),'-dpng')


%% chl
fig=figure; hold on;
i = 1:4
temp_data = nanmean(chl_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(chl_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('1to4','6to9')
print(fig,strcat('chl_b_sp_mean.png'),'-dpng')

fig=figure; hold on;
i = 2:4
temp_data = nanmean(chl_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data,'linew',2);
% plot(2:3:264,nanmean(yoonja.chl_sur(2:6,2:3:264),1),'color','b','linew',2);
xticklabels(1997:2019);
xticks(1:23);
xlim([0 23]); grid on;
ylim([0 14]);
ylabel('chl (ug/L)')
title('chl bottom')
xtickangle(45);
i = 6:9
temp_data = nanmean(chl_b_yr(i,:),1);
no_nan_data= temp_data(~isnan(temp_data));
no_nan_ind=find(isnan(temp_data)==0);
plot(no_nan_ind,no_nan_data ,'linew',2);
legend('2to4','6to9')
print(fig,strcat('chl_b_sp_mean(minus1).png'),'-dpng')


%%
fig=figure; hold on;
% plot((gy_csh_in_re)','color',[.5,.5,.5]);
plot(2:3:264,regime_temp(sp_gy(I(j)),2:3:264),'color','r','linew',2);
plot(2:3:264,yoonja.temp_sur(j,2:3:264),'color','b','linew',2);
xticklabels(1997:2018);
xticks(1:12:264);
xlim([0 264]); grid on;
ylim([0 30])
ylabel('temp (^oC)')
title('temp surface')
xtickangle(45);
legend('cshwa','Pf.Kang','fontsize',15);
print(fig,strcat('label.png'),'-dpng')