close all; clear; clc; 
[raw1,txt1]=xlsread('yeosu_1980.xlsx','sheet2','');
% for i= 1:8784 %366*24
%     elev(i)=nanmean(raw1(1+(i-1)*60:60*i,1),1);
% end   
elev = raw1(:,6)./100; 
save('Yeosu_tide_2001.mat','elev')

close all; clear; clc;
load Yeosu_tide_2001.mat
z=ncread('merg_z_hr.nc', 'zeta');

z_st = squeeze(z(74,24,:));

figure; plot(z_st); hold on;
plot(elev(1:length(z_st))-mean(elev),'r')
ylim([-1.2 1.2])
ylabel('elevation(m)','fontsize',12,'fontweight','bold')
xlabel('hours','fontsize',12,'fontweight','bold'); grid on;
set(gca,'fontsize',12,'fontweight','bold');

% [raw1,txt1]=xlsread('yeosu_1980.xlsx','sheet1','');
% % for i= 1:8784 %366*24
% %     elev(i)=nanmean(raw1(1+(i-1)*60:60*i,1),1);
% % end   
% elev = raw1(:,6); 
% save('Yeosu_tide_1980.mat','elev');

close all; clear; clc; 
%% find nearest point on Ye
grd_f = 'grid_gy_v11_s.nc';
lon = ncread(grd_f, 'lon_rho');
lat = ncread(grd_f, 'lat_rho');
diff=[lon(1:end,1)] - [127.765556];
min(abs(diff))
find(min(abs(diff))==diff)
diff=[lat(1,1:end)] - [34.74721840];
min(abs(diff))
find(min(abs(diff))==diff)

j=0;
% NAN(252,176,8760)
for i = 1:8760 %% for 999 file
% for i = 1000:8760
        filename = strcat('ocean_his_',num2str(i,'%04i'),'*.nc'); %% for 999 file
%         filename = strcat('pohang_avg_',num2str(i),'*.nc'); %% from 1000 ~ end file
        d_filename = dir(filename);       
        j = j+1;
        for k = 1:length(d_filename)
            temp = ncread(d_filename(k).name, 'zeta');
              j
            total_temp(:,j) = squeeze(temp(74,24));
        end         
end

% total = permute(total_temp,[2 1 3 4]);
% temp = total;
% save('temp_pohang_gg1','temp','-v7.3'); %759 loop >> 999 file
% save('temp_pohang_2','temp','-v7.3'); % 1881 loop >> 2880 filer



% ����Ҹ�?����
% ����:34.74721840
% �浵:127.765556
load Yeosu_tide_1980.mat

ssh = total_temp;
ano_elev=elev - mean(elev);

% for i = 1:size(total_temp,3)
% ssh(i)=griddata(lon,lat,squeeze(total_temp(:,:,i)),127.765556,34.74721840,'nearest');
% i
% end
% save('elev_to_269hr.mat');
save('elev_to_781hr.mat');



figure
plot(9:322,ssh(1:end))
hold on
plot((elev(1:314)-nanmean(elev))./100,'r')

figure
plot(ssh-mean(ssh))
hold on
% plot((elev(1:1169)-nanmean(elev))./100,'r')
plot((elev(7*24-2:end)-nanmean(elev))./100,'r')

figure
plot(ssh(8*24:end))
hold on
% plot((elev(1:1169)-nanmean(elev))./100,'r')
plot((elev(1:1169)-nanmean(elev(1:1169)))./100,'r')

figure
plot(ssh(1:end))
hold on
% plot((elev(1:1169)-nanmean(elev))./100,'r')
plot((elev(1:1169)-nanmean(elev(1:1169)))./100,'r')



figure
plot(ssh(1:end),'linew',2)
hold on
% plot((elev(1:1169)-nanmean(elev))./100,'r')
plot((ano_elev(1:742))./100,'r','linew',1.5)
xlabel(gca,'elevation(m)'); ylabel('time(hr)');
legend(gca,'model', 'TG'); grid on;
set(gca,'fontsize',18, 'fontweight','bold');
