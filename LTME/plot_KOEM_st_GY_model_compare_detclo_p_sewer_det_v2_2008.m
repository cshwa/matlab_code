close all; clear; clc;

cd F:\ROMS\roms_tools\Run\
start

cd D:\Àå±â»ıÅÂ\Dynamic\KOEM
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
mod_file='D:\Àå±â»ıÅÂ\Dynamic\result\2008\2008_mp_p_sewer_det_f.nc';
mod_no3=ncread(mod_file, 'NO3');
mod_nh4=ncread(mod_file, 'NH4');
mod_chl=ncread(mod_file, 'chlorophyll');
mod_do=ncread(mod_file, 'oxygen');
mod_temp=ncread(mod_file, 'temp');
mod_salt=ncread(mod_file, 'salt');
mod_po4=ncread(mod_file, 'tPO4');
mod_zoo=ncread(mod_file, 'zooplankton');
mod_phy=ncread(mod_file, 'phytoplankton');

% missing value will be NaN;
mod_no3(mod_no3 > 10^25) =NaN; 
mod_nh4(mod_nh4 > 10^25) =NaN; 
mod_do(mod_do > 10^25) =NaN; 
mod_chl(mod_chl > 10^25) =NaN; 
mod_temp(mod_temp > 10^25) =NaN; 
mod_salt(mod_salt > 10^25) =NaN; 
mod_po4(mod_po4 > 10^25) =NaN; 
mod_zoo(mod_zoo > 10^25) =NaN; 
mod_phy(mod_phy > 10^25) =NaN; 

% if leapyear(2002) ==  0
%     mod_no3(:,:,:,366) =[]; 
%     mod_nh4(:,:,:,366) =[]; 
%     mod_do(:,:,:,366) =[]; 
%     mod_chl(:,:,:,366) =[];  
%     mod_temp(:,:,:,366) =[];  
%     mod_salt(:,:,:,366) =[];  
%     mod_po4(:,:,:,366) =[]; 
%     mod_zoo(:,:,:,366) =[];  
%     mod_phy(:,:,:,366) =[];  
% end

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
for i = 2008:2008
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
clearvars mod_m_*
if size(mod_no3,4) >= 365
    for i =1:12
    if i ==1 
        mod_m_no3(:,:,:,i)=mean(mod_no3(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_nh4(:,:,:,i)=mean(mod_nh4(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_do(:,:,:,i)=mean(mod_do(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_chl(:,:,:,i)=mean(mod_chl(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_temp(:,:,:,i)=mean(mod_temp(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_salt(:,:,:,i)=mean(mod_salt(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_po4(:,:,:,i)=mean(mod_po4(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_zoo(:,:,:,i)=mean(mod_zoo(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_phy(:,:,:,i)=mean(mod_phy(:,:,:,1:eom_d_each(1,i)),4);
    else
        mod_m_no3(:,:,:,i)=mean(mod_no3(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_nh4(:,:,:,i)=mean(mod_nh4(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_do(:,:,:,i)=mean(mod_do(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_chl(:,:,:,i)=mean(mod_chl(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_temp(:,:,:,i)=mean(mod_temp(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_salt(:,:,:,i)=mean(mod_salt(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_po4(:,:,:,i)=mean(mod_po4(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_zoo(:,:,:,i)=mean(mod_zoo(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_phy(:,:,:,i)=mean(mod_phy(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
    end
    end
else
    for i =1:10
    if i ==1 
        mod_m_no3(:,:,:,i)=mean(mod_no3(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_nh4(:,:,:,i)=mean(mod_nh4(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_do(:,:,:,i)=mean(mod_do(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_chl(:,:,:,i)=mean(mod_chl(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_temp(:,:,:,i)=mean(mod_temp(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_salt(:,:,:,i)=mean(mod_salt(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_po4(:,:,:,i)=mean(mod_po4(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_zoo(:,:,:,i)=mean(mod_zoo(:,:,:,1:eom_d_each(1,i)),4);
        mod_m_phy(:,:,:,i)=mean(mod_phy(:,:,:,1:eom_d_each(1,i)),4);
    else
        mod_m_no3(:,:,:,i)=mean(mod_no3(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_nh4(:,:,:,i)=mean(mod_nh4(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_do(:,:,:,i)=mean(mod_do(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_chl(:,:,:,i)=mean(mod_chl(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_temp(:,:,:,i)=mean(mod_temp(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_salt(:,:,:,i)=mean(mod_salt(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_po4(:,:,:,i)=mean(mod_po4(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_zoo(:,:,:,i)=mean(mod_zoo(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
        mod_m_phy(:,:,:,i)=mean(mod_phy(:,:,:,eom_d_each(1,i-1)+1:eom_d_each(1,i)),4);
    end
    end
end
    
% obs_tdx = 49:60; %2001 %year
% obs_tdx = 97:108; %2005 %year
% obs_tdx = 109:120; %2006 %year
% obs_tdx = 121:132; %2007 %year
obs_tdx = 133:144; %2008 %year
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

% compare plot KOEM vs. MODEL
n_2d_x = near_point_2d(:,1).*(ind_out'./ind_out'); n_2d_y = near_point_2d(:,2).*(ind_out'./ind_out');
n_2d_x(isnan(n_2d_x))=[]; n_2d_y(isnan(n_2d_y))=[];

%%%%%% limit check %%%%%
%%%%%% surf %%%%%
max(max(max(squeeze(mod_no3(n_2d_x,n_2d_y,20,:)))))
min(min(min(squeeze(mod_no3(n_2d_x,n_2d_y,20,:)))))

max(max(max(squeeze(mod_nh4(n_2d_x,n_2d_y,20,:)))))
min(min(min(squeeze(mod_nh4(n_2d_x,n_2d_y,20,:)))))

max(max(max(squeeze(mod_do(n_2d_x,n_2d_y,20,:)))))
min(min(min(squeeze(mod_do(n_2d_x,n_2d_y,20,:)))))

max(max(max(squeeze(mod_temp(n_2d_x,n_2d_y,20,:)))))
min(min(min(squeeze(mod_temp(n_2d_x,n_2d_y,20,:)))))

max(max(max(squeeze(mod_salt(n_2d_x,n_2d_y,20,:)))))
min(min(min(squeeze(mod_salt(n_2d_x,n_2d_y,20,:)))))

max(max(max(squeeze(mod_chl(n_2d_x,n_2d_y,20,:)))))
min(min(min(squeeze(mod_chl(n_2d_x,n_2d_y,20,:)))))

%%%%%% bot %%%%%
max(max(max(squeeze(mod_no3(n_2d_x,n_2d_y,1,:)))))
min(min(min(squeeze(mod_no3(n_2d_x,n_2d_y,1,:)))))

max(max(max(squeeze(mod_nh4(n_2d_x,n_2d_y,1,:)))))
min(min(min(squeeze(mod_nh4(n_2d_x,n_2d_y,1,:)))))

max(max(max(squeeze(mod_do(n_2d_x,n_2d_y,1,:)))))
min(min(min(squeeze(mod_do(n_2d_x,n_2d_y,1,:)))))

max(max(max(squeeze(mod_temp(n_2d_x,n_2d_y,1,:)))))
min(min(min(squeeze(mod_temp(n_2d_x,n_2d_y,1,:)))))

max(max(max(squeeze(mod_salt(n_2d_x,n_2d_y,1,:)))))
min(min(min(squeeze(mod_salt(n_2d_x,n_2d_y,1,:)))))

max(max(max(squeeze(mod_chl(n_2d_x,n_2d_y,1,:)))))
min(min(min(squeeze(mod_chl(n_2d_x,n_2d_y,1,:)))))

% %%%%%%%%%%
% return
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_no3(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_no3(i,:)./14,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NO3']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%         close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_nh4(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_nh4(i,:)./14,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NH4']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_NH4_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off; close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_do(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_do(i,:).*0.7.*44.661,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL DO']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('DO (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_DO_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_temp(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_temp(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL temp']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('temp (^oC)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_temp_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_salt(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_salt(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL salt']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_salt_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_chl(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b','linew',2);
%         plot(obs_chl(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL chla']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_chl_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%         close;
%    end
% end
% 
% %% bottom daily
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_no3(near_point_2d(i,1),near_point_2d(i,2),1,:)),'b','linew',2);
%         plot(obs_no3_b(i,:)./14,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NO3']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('NO3 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_b_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%         close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_nh4(near_point_2d(i,1),near_point_2d(i,2),1,:)),'b','linew',2);
%         plot(obs_nh4_b(i,:)./14,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NH4']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('NH4 (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_b_NH4_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off; close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_do(near_point_2d(i,1),near_point_2d(i,2),1,:)),'b','linew',2);
%         plot(obs_do_b(i,:).*0.7.*44.661,'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL DO']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('DO (umol/m^3)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_DO_KOEM_b_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_temp(near_point_2d(i,1),near_point_2d(i,2),1,:)),'b','linew',2);
%         plot(obs_temp_b(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL temp']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('temp (^oC)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_temp_b_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
% %         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_salt(near_point_2d(i,1),near_point_2d(i,2),1,:)),'b','linew',2);
%         plot(obs_salt_b(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL salt']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('salt (psu)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_salt_b_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;close;
%    end
% end
% 
% %plot daily
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
%         plot(squeeze(mod_chl(near_point_2d(i,1),near_point_2d(i,2),1,:)),'b','linew',2);
%         plot(obs_chl_b(i,:),'r-','linew',2); xlim([1 365]);
%         title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL chla']);
%         xlabel('time(days on 2002)','fontsize',13)
%         ylabel('chl (ug/L)','fontsize',13)
%         grid on
%         set(gca,'fontsize',13)
%         name_gf= name_tag{i};
%         print(fig,strcat('2001_daily_chl_b_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%         close;
%    end
% end
% 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % plot monthly
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        figure; hold on;
%         plot(squeeze(mod_m_no3(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b');
%         plot(regime_no3(i,49:60)./14,'ro'); xlim([1 12])
%         title([name_tag{i},'- 2001 monthly KOEM OBS vs. MODEL']);
%         print(fig,strcat('2001_mon_no3_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%    end
% end
% 
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        figure; hold on;
%         plot(squeeze(mod_m_nh4(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b');
%         plot(regime_nh4(i,49:60)./14,'ro'); xlim([1 12])
%         title([name_tag{i},'- 2001 monthly KOEM OBS vs. MODEL']);
%         print(fig,strcat('2001_mon_nh4_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%    end
% end
% 
% 
% % plot monthly
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        figure; hold on;
%         plot(squeeze(mod_m_chl(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b');
%         plot(regime_no3(i,49:60)./14,'ro'); xlim([1 12])
%         title([name_tag{i},'- 2001 monthly KOEM OBS vs. MODEL']);
%         print(fig,strcat('2001_mon_temp_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%    end
% end
% 
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        figure; hold on;
%         plot(squeeze(mod_m_no3(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b');
%         plot(regime_no3(i,49:60)./14,'ro'); xlim([1 12])
%         title([name_tag{i},'- 2001 monthly KOEM OBS vs. MODEL']);
%         print(fig,strcat('2001_mon_salt_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%    end
% end
% 
% 
% for i = 1:65
% %   for   i = 4:4
%    if ind_out(i) == 1
%        figure; hold on;
%         plot(squeeze(mod_m_chl(near_point_2d(i,1),near_point_2d(i,2),20,:)),'b');
%         plot(regime_chl(i,49:60),'ro'); xlim([1 12])
%         title([name_tag{i},'- 2001 monthly KOEM OBS vs. MODEL chla']);
%         print(fig,strcat('2001_mon_chla_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
%         hold off;
%    end
% end
% 
% return
% 
% 
% sp_mean_no3_yr
% sp_mean_temp_yr
% sp_mean_salt_yr
% %% salt
% clearvars yp_w_salt
% color_pic = lines(size(regime_salt_yr,1));
% marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<',...
%     'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<'};
% xp = 1:22;
% j=0
% figure; hold on;
% for i = 1:size(sp_mean_salt_yr,1)
%   clearvars reg_data_salt xp_w_salt pf_w_salt
% reg_data_salt = sp_mean_salt_yr(i,:);
% if isnan(reg_data_salt(1)) == 0 
%     j = j+1;
%     xp_w_salt = find(isnan(reg_data_salt)==0);
%     pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
%     yp_w_salt(i,:) = polyval(pf_w_salt,xp);
%     scatter(1:22,sp_mean_salt_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
%     plot(1:22, yp_w_salt(i,:),'color',color_pic(i,:));
%     sp_coeff_salt(i,:) = pf_w_salt;
%     sp_temp_case(j) = i;
% end
% hold on
% end
% xlabel('time(year)','fontsize',13)
% ylabel('salt (mg/m^3)','fontsize',13)
% set(gca,'xtick',[1:2:22]);
% set(gca,'xlim',[1 22]);
% set(gca,'xticklabel',1997:2:2018);
% title('KOEM-Ç¥Ãş¿°ºĞ ¿¬Æò±Õ(°ø°£Æò±Õ)','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% % ylim([32 35])
% 
% 
% %% no3
% clearvars yp_w_no3
% j=0
% figure; hold on;
% for i = 1:size(sp_mean_no3_yr,1)
%   clearvars reg_data_no3 xp_w_no3 pf_w_no3
% reg_data_no3 = sp_mean_no3_yr(i,:);
% if isnan(reg_data_no3(1)) == 0 
%     j = j+1;
%     xp_w_no3 = find(isnan(reg_data_no3)==0);
%     pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
%     yp_w_no3(i,:) = polyval(pf_w_no3,xp);
%     scatter(1:22,sp_mean_no3_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
%     plot(1:22, yp_w_no3(i,:),'color',color_pic(i,:));
%     sp_coeff_no3(i,:) = pf_w_no3;
%     sp_temp_case(j) = i;
% end
% hold on
% end
% xlabel('time(year)','fontsize',13)
% ylabel('NO3 (ug/L)','fontsize',13)
% set(gca,'xtick',[1:2:22]);
% set(gca,'xlim',[1 22]);
% set(gca,'xticklabel',1997:2:2018);
% title('KOEM°üÃø-Ç¥ÃşÁú»ê¿° ¿¬Æò±Õ(°ø°£Æò±Õ)','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% % ylim([32 35])
% % legend('205-01','205-02','205-03','205-04','205-05')
% 
% %temp
% clearvars yp_w_temp
% j=0
% figure; hold on;
% for i = 1:size(sp_mean_temp_yr,1)
%   clearvars reg_data_temp xp_w_temp pf_w_temp
% reg_data_temp = sp_mean_temp_yr(i,:);
% if isnan(reg_data_temp(1)) == 0 
%     j = j+1;
%     xp_w_temp = find(isnan(reg_data_temp)==0);
%     pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
%     yp_w_temp(i,:) = polyval(pf_w_temp,xp);
%     scatter(1:22,sp_mean_temp_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
%     plot(1:22, yp_w_temp(i,:),'color',color_pic(i,:));
%     sp_coeff_temp(i,:) = pf_w_temp;
%     sp_temp_case(j) = i;
% end
% hold on
% end
% xlabel('time(year)','fontsize',13)
% ylabel('temperature (^oC)','fontsize',13)
% set(gca,'xtick',[1:2:22]);
% set(gca,'xlim',[1 22]);
% set(gca,'xticklabel',1997:2:2018);
% title('KOEM°üÃø-Ç¥Ãş¼ö¿Â ¿¬Æò±Õ','fontsize',13)
% grid on
% set(gca,'fontsize',13)
% ylim([13 20])

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

squeeze(nanmean(nanmean(mod_no3(near_point_2d(sp_gy,1), near_point_2d(sp_gy,2),20,:),1),2));
clearvars sp_std_*
for j = 1:length(sp_gy)
    sp_gy_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_no3_mid(j,:)=squeeze(mod_no3(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),10,:));
    sp_gy_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_nh4_mid(j,:)=squeeze(mod_nh4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),10,:));
    sp_gy_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_chl_mid(j,:)=squeeze(mod_chl(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),10,:));
    sp_gy_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_temp_mid(j,:)=squeeze(mod_temp(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),10,:));
    sp_gy_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_salt_mid(j,:)=squeeze(mod_salt(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),10,:));
    sp_gy_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_do(j,:)=squeeze(mod_do(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_do_mid(j,:)=squeeze(mod_do(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),10,:)); 
    sp_gy_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_po4(j,:)=squeeze(mod_po4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_po4_mid(j,:)=squeeze(mod_po4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),10,:)); 
    sp_gy_po4_b(j,:)=squeeze(mod_po4(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));   
    sp_gy_zoo(j,:)=squeeze(mod_zoo(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_zoo_b(j,:)=squeeze(mod_zoo(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
    sp_gy_phy(j,:)=squeeze(mod_phy(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),20,:));
    sp_gy_phy_b(j,:)=squeeze(mod_phy(near_point_2d(sp_gy(j),1), near_point_2d(sp_gy(j),2),1,:));
end

for j = 1:length(sp_s_gy)
    sp_sgy_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_no3_mid(j,:)=squeeze(mod_no3(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),10,:));
    sp_sgy_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_nh4_mid(j,:)=squeeze(mod_nh4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),10,:));
    sp_sgy_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_chl_mid(j,:)=squeeze(mod_chl(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),10,:));
    sp_sgy_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_temp_mid(j,:)=squeeze(mod_temp(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),10,:));
    sp_sgy_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_salt_mid(j,:)=squeeze(mod_salt(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),10,:));
    sp_sgy_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_do(j,:)=squeeze(mod_do(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_do_mid(j,:)=squeeze(mod_do(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),10,:)); 
    sp_sgy_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_po4(j,:)=squeeze(mod_po4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_po4_mid(j,:)=squeeze(mod_po4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),10,:)); 
    sp_sgy_po4_b(j,:)=squeeze(mod_po4(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_zoo(j,:)=squeeze(mod_zoo(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_zoo_b(j,:)=squeeze(mod_zoo(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
    sp_sgy_phy(j,:)=squeeze(mod_phy(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),20,:));
    sp_sgy_phy_b(j,:)=squeeze(mod_phy(near_point_2d(sp_s_gy(j),1), near_point_2d(sp_s_gy(j),2),1,:));
end

for j = 1:length(sp_e_gy)
    sp_egy_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_no3_mid(j,:)=squeeze(mod_no3(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),10,:));
    sp_egy_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_nh4_mid(j,:)=squeeze(mod_nh4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),10,:));
    sp_egy_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_chl_mid(j,:)=squeeze(mod_chl(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),10,:));
    sp_egy_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_temp_mid(j,:)=squeeze(mod_temp(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),10,:));
    sp_egy_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_salt_mid(j,:)=squeeze(mod_salt(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),10,:));
    sp_egy_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_do(j,:)=squeeze(mod_do(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_do_mid(j,:)=squeeze(mod_do(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),10,:));
    sp_egy_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));  
    sp_egy_po4(j,:)=squeeze(mod_po4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_po4_mid(j,:)=squeeze(mod_po4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),10,:));
    sp_egy_po4_b(j,:)=squeeze(mod_po4(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_zoo(j,:)=squeeze(mod_zoo(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_zoo_b(j,:)=squeeze(mod_zoo(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
    sp_egy_phy(j,:)=squeeze(mod_phy(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),20,:));
    sp_egy_phy_b(j,:)=squeeze(mod_phy(near_point_2d(sp_e_gy(j),1), near_point_2d(sp_e_gy(j),2),1,:));
end

for j = 1:length(sp_jj)
    sp_jj_no3(j,:)=squeeze(mod_no3(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_no3_mid(j,:)=squeeze(mod_no3(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),10,:));
    sp_jj_no3_b(j,:)=squeeze(mod_no3(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_nh4(j,:)=squeeze(mod_nh4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_nh4_mid(j,:)=squeeze(mod_nh4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),10,:));
    sp_jj_nh4_b(j,:)=squeeze(mod_nh4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_chl(j,:)=squeeze(mod_chl(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_chl_mid(j,:)=squeeze(mod_chl(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),10,:));
    sp_jj_chl_b(j,:)=squeeze(mod_chl(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_temp(j,:)=squeeze(mod_temp(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_temp_mid(j,:)=squeeze(mod_temp(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),10,:));
    sp_jj_temp_b(j,:)=squeeze(mod_temp(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_salt(j,:)=squeeze(mod_salt(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_salt_mid(j,:)=squeeze(mod_salt(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),10,:));
    sp_jj_salt_b(j,:)=squeeze(mod_salt(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_do(j,:)=squeeze(mod_do(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_do_mid(j,:)=squeeze(mod_do(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),10,:));
    sp_jj_do_b(j,:)=squeeze(mod_do(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));  
    sp_jj_po4(j,:)=squeeze(mod_po4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_po4_mid(j,:)=squeeze(mod_po4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),10,:));
    sp_jj_po4_b(j,:)=squeeze(mod_po4(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_zoo(j,:)=squeeze(mod_zoo(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_zoo_b(j,:)=squeeze(mod_zoo(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    sp_jj_phy(j,:)=squeeze(mod_phy(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),20,:));
    sp_jj_phy_b(j,:)=squeeze(mod_phy(near_point_2d(sp_jj(j),1), near_point_2d(sp_jj(j),2),1,:));
    
end


for i = 1:size(sp_gy_no3,2)
    %std
    std_gy_no3(i)=nanstd(squeeze(sp_gy_no3(:,i)));
    std_gy_no3_mid(i)=nanstd(squeeze(sp_gy_no3_mid(:,i)));
    std_gy_no3_b(i)=nanstd(squeeze(sp_gy_no3_b(:,i)));
    std_gy_nh4(i)=nanstd(squeeze(sp_gy_nh4(:,i)));
    std_gy_nh4_mid(i)=nanstd(squeeze(sp_gy_nh4_mid(:,i)));
    std_gy_nh4_b(i)=nanstd(squeeze(sp_gy_nh4_b(:,i)));
    std_gy_chl(i)=nanstd(squeeze(sp_gy_chl(:,i)));
    std_gy_chl_mid(i)=nanstd(squeeze(sp_gy_chl_mid(:,i)));
    std_gy_chl_b(i)=nanstd(squeeze(sp_gy_chl_b(:,i)));
    std_gy_temp(i)=nanstd(squeeze(sp_gy_temp(:,i)));
    std_gy_temp_mid(i)=nanstd(squeeze(sp_gy_temp_mid(:,i)));
    std_gy_temp_b(i)=nanstd(squeeze(sp_gy_temp_b(:,i)));
    std_gy_salt(i)=nanstd(squeeze(sp_gy_salt(:,i)));
    std_gy_salt_mid(i)=nanstd(squeeze(sp_gy_salt_mid(:,i)));
    std_gy_salt_b(i)=nanstd(squeeze(sp_gy_salt_b(:,i)));
    std_gy_do(i)=nanstd(squeeze(sp_gy_do(:,i)));
    std_gy_do_mid(i)=nanstd(squeeze(sp_gy_do_mid(:,i)));
    std_gy_do_b(i)=nanstd(squeeze(sp_gy_do_b(:,i)));
    
    std_gy_zoo(i)=nanstd(squeeze(sp_gy_zoo(:,i)));
    std_gy_zoo_b(i)=nanstd(squeeze(sp_gy_zoo_b(:,i)));
    std_gy_phy(i)=nanstd(squeeze(sp_gy_phy(:,i)));
    std_gy_phy_b(i)=nanstd(squeeze(sp_gy_phy_b(:,i)));


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
    std_sgy_zoo(i)=nanstd(squeeze(sp_sgy_zoo(:,i)));
    std_sgy_zoo_b(i)=nanstd(squeeze(sp_sgy_zoo_b(:,i)));
    std_sgy_phy(i)=nanstd(squeeze(sp_sgy_phy(:,i)));
    std_sgy_phy_b(i)=nanstd(squeeze(sp_sgy_phy_b(:,i)));
    

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
    std_egy_zoo(i)=nanstd(squeeze(sp_egy_zoo(:,i)));
    std_egy_zoo_b(i)=nanstd(squeeze(sp_egy_zoo_b(:,i)));
    std_egy_phy(i)=nanstd(squeeze(sp_egy_phy(:,i)));
    std_egy_phy_b(i)=nanstd(squeeze(sp_egy_phy_b(:,i)));
    
    std_egy_do(i)=nanstd(squeeze(sp_egy_do(:,i)));
    std_egy_do_b(i)=nanstd(squeeze(sp_egy_do_b(:,i)));
    std_egy_no3_mid(i)=nanstd(squeeze(sp_egy_no3_mid(:,i)));
    std_egy_nh4_mid(i)=nanstd(squeeze(sp_egy_nh4_mid(:,i)));
    std_egy_chl_mid(i)=nanstd(squeeze(sp_egy_chl_mid(:,i)));
    std_egy_temp_mid(i)=nanstd(squeeze(sp_egy_temp_mid(:,i)));
    std_egy_salt_mid(i)=nanstd(squeeze(sp_egy_salt_mid(:,i)));
    std_egy_do_mid(i)=nanstd(squeeze(sp_egy_do_mid(:,i)));

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
    std_jj_zoo(i)=nanstd(squeeze(sp_jj_zoo(:,i)));
    std_jj_zoo_b(i)=nanstd(squeeze(sp_jj_zoo_b(:,i)));
    std_jj_phy(i)=nanstd(squeeze(sp_jj_phy(:,i)));
    std_jj_phy_b(i)=nanstd(squeeze(sp_jj_phy_b(:,i)));
    
    std_jj_no3_mid(i)=nanstd(squeeze(sp_jj_no3_mid(:,i)));
    std_jj_nh4_mid(i)=nanstd(squeeze(sp_jj_nh4_mid(:,i)));
    std_jj_chl_mid(i)=nanstd(squeeze(sp_jj_chl_mid(:,i)));
    std_jj_temp_mid(i)=nanstd(squeeze(sp_jj_temp_mid(:,i)));
    std_jj_salt_mid(i)=nanstd(squeeze(sp_jj_salt_mid(:,i)));
    std_jj_do_mid(i)=nanstd(squeeze(sp_jj_do_mid(:,i)));
    
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
    gy_zoo(i)=nanmean(squeeze(sp_gy_zoo(:,i)));
    gy_zoo_b(i)=nanmean(squeeze(sp_gy_zoo_b(:,i)));
    gy_phy(i)=nanmean(squeeze(sp_gy_phy(:,i)));
    gy_phy_b(i)=nanmean(squeeze(sp_gy_phy_b(:,i)));

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
    sgy_zoo(i)=nanmean(squeeze(sp_sgy_zoo(:,i)));
    sgy_zoo_b(i)=nanmean(squeeze(sp_sgy_zoo_b(:,i)));
    sgy_phy(i)=nanmean(squeeze(sp_sgy_phy(:,i)));
    sgy_phy_b(i)=nanmean(squeeze(sp_sgy_phy_b(:,i)));

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
    egy_zoo(i)=nanmean(squeeze(sp_egy_zoo(:,i)));
    egy_zoo_b(i)=nanmean(squeeze(sp_egy_zoo_b(:,i)));
    egy_phy(i)=nanmean(squeeze(sp_egy_phy(:,i)));
    egy_phy_b(i)=nanmean(squeeze(sp_egy_phy_b(:,i)));

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
    jj_zoo(i)=nanmean(squeeze(sp_jj_zoo(:,i)));
    jj_zoo_b(i)=nanmean(squeeze(sp_jj_zoo_b(:,i)));
    jj_phy(i)=nanmean(squeeze(sp_jj_phy(:,i)));
    jj_phy_b(i)=nanmean(squeeze(sp_jj_phy_b(:,i)));
    
    gy_no3_mid(i)=nanmean(squeeze(sp_gy_no3_mid(:,i)));
    gy_nh4_mid(i)=nanmean(squeeze(sp_gy_nh4_mid(:,i)));
    gy_chl_mid(i)=nanmean(squeeze(sp_gy_chl_mid(:,i)));
    gy_temp_mid(i)=nanmean(squeeze(sp_gy_temp_mid(:,i)));
    gy_salt_mid(i)=nanmean(squeeze(sp_gy_salt_mid(:,i)));
    gy_do_mid(i)=nanmean(squeeze(sp_gy_do_mid(:,i)));
    
    
    sgy_no3_mid(i)=nanmean(squeeze(sp_sgy_no3_mid(:,i)));
    sgy_nh4_mid(i)=nanmean(squeeze(sp_sgy_nh4_mid(:,i)));
    sgy_chl_mid(i)=nanmean(squeeze(sp_sgy_chl_mid(:,i)));
    sgy_temp_mid(i)=nanmean(squeeze(sp_sgy_temp_mid(:,i)));
    sgy_salt_mid(i)=nanmean(squeeze(sp_sgy_salt_mid(:,i)));
    sgy_do_mid(i)=nanmean(squeeze(sp_sgy_do_mid(:,i)));

   
    egy_no3_mid(i)=nanmean(squeeze(sp_egy_no3_mid(:,i)));
    egy_nh4_mid(i)=nanmean(squeeze(sp_egy_nh4_mid(:,i)));
    egy_chl_mid(i)=nanmean(squeeze(sp_egy_chl_mid(:,i)));
    egy_temp_mid(i)=nanmean(squeeze(sp_egy_temp_mid(:,i)));
    egy_salt_mid(i)=nanmean(squeeze(sp_egy_salt_mid(:,i)));
    egy_do_mid(i)=nanmean(squeeze(sp_egy_do_mid(:,i)));

    
    jj_no3_mid(i)=nanmean(squeeze(sp_jj_no3_mid(:,i)));
    jj_nh4_mid(i)=nanmean(squeeze(sp_jj_nh4_mid(:,i)));
    jj_chl_mid(i)=nanmean(squeeze(sp_jj_chl_mid(:,i)));
    jj_temp_mid(i)=nanmean(squeeze(sp_jj_temp_mid(:,i)));
    jj_salt_mid(i)=nanmean(squeeze(sp_jj_salt_mid(:,i)));
    jj_do_mid(i)=nanmean(squeeze(sp_jj_do_mid(:,i)));
    
    %po4
    std_gy_po4(i)=nanstd(squeeze(sp_gy_po4(:,i)));
    std_gy_po4_mid(i)=nanstd(squeeze(sp_gy_po4_mid(:,i)));
    std_gy_po4_b(i)=nanstd(squeeze(sp_gy_po4_b(:,i)));
    
    std_sgy_po4(i)=nanstd(squeeze(sp_sgy_po4(:,i)));
    std_sgy_po4_b(i)=nanstd(squeeze(sp_sgy_po4_b(:,i)));
    
    std_egy_po4(i)=nanstd(squeeze(sp_egy_po4(:,i)));
    std_egy_po4_b(i)=nanstd(squeeze(sp_egy_po4_b(:,i)));
    std_egy_po4_mid(i)=nanstd(squeeze(sp_egy_po4_mid(:,i)));
    
    
    std_jj_po4(i)=nanstd(squeeze(sp_jj_po4(:,i)));
    std_jj_po4_b(i)=nanstd(squeeze(sp_jj_po4_b(:,i)));
    
    
     gy_po4(i)=nanmean(squeeze(sp_gy_po4(:,i)));
    gy_po4_b(i)=nanmean(squeeze(sp_gy_po4_b(:,i)));


     sgy_po4(i)=nanmean(squeeze(sp_sgy_po4(:,i)));
    sgy_po4_b(i)=nanmean(squeeze(sp_sgy_po4_b(:,i)));
    
    egy_po4(i)=nanmean(squeeze(sp_egy_po4(:,i)));
    egy_po4_b(i)=nanmean(squeeze(sp_egy_po4_b(:,i)));
    
        jj_po4(i)=nanmean(squeeze(sp_jj_po4(:,i)));
    jj_po4_b(i)=nanmean(squeeze(sp_jj_po4_b(:,i)));
    
    gy_po4_mid(i)=nanmean(squeeze(sp_gy_po4_mid(:,i)));
    sgy_po4_mid(i)=nanmean(squeeze(sp_sgy_po4_mid(:,i)));
    egy_po4_mid(i)=nanmean(squeeze(sp_egy_po4_mid(:,i)));
     jj_po4_mid(i)=nanmean(squeeze(sp_jj_po4_mid(:,i)));
  
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

cd D:\Àå±â»ıÅÂ\Dynamic\result\2008
% save('2002_mp_p_sewer_det_result.mat','jj_*', 'gy_*', 'sgy_*', 'egy_*','obm*','obs*');     
save('2008_mp_p_sewer_det_f_result.mat','-v7.3');     

return
%% GY plot
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_zoo,'b','linew',2);
        plot(gy_phy,'r','linew',2);
        plot(gy_zoo_b,'g--','linew',2);
        plot(gy_phy_b,'m--','linew',2);
        
        legend('zoo','phy','zoo_b_o_t','phy_b_o_t')
        
%         plot(obm_gy_po4./14,'r-','linew',2); 
        xlim([1 365]);
%         plot(1:length(obm_gy_po4),obm_gy_po4./14 + obs_std_gy_po4./14,'m-','linew',2);
%         plot(1:length(obm_gy_po4),obm_gy_po4./14 - obs_std_gy_po4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL Zoo & Phy']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('concentration (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)

  %% temperature      
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);
        plot(gy_temp_b,'r','linew',2);
        nonan_data_plt = find(isnan(obm_gy_temp)==0);
        nonan_data_plt_b = find(isnan(obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        
        nonan_diff_b=nonan_data_plt_b(2:end) - nonan_data_plt_b(1:end-1);
        discon_p_b = find(nonan_diff_b ~= 1); % find discontinuous point
        discon_p_b(4) = length(nonan_data_plt_b); %end point
        
        upper_bound_plt = obm_gy_temp + obs_std_gy_temp;
        lower_bound_plt = obm_gy_temp - obs_std_gy_temp;
        upper_bound_plt_b = obm_gy_temp_b + obs_std_gy_temp_b;
        lower_bound_plt_b = obm_gy_temp_b - obs_std_gy_temp_b;
      
        plot(obm_gy_temp,'b-','linew',2); xlim([1 365]);
        plot(obm_gy_temp_b,'r-','linew',2);
        plot(1:length(obm_gy_temp),upper_bound_plt,'b-','linew',2);
        plot(1:length(obm_gy_temp),lower_bound_plt,'b-','linew',2);
        plot(1:length(obm_gy_temp_b),upper_bound_plt_b,'r-','linew',2);
        plot(1:length(obm_gy_temp_b),lower_bound_plt_b,'r-','linew',2);
        
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'b');
        patch([nonan_data_plt_b(1:discon_p_b(1)) fliplr(nonan_data_plt_b(1:discon_p_b(1)))], [upper_bound_plt_b(nonan_data_plt_b(1:discon_p_b(1))) lower_bound_plt_b(nonan_data_plt_b(1:discon_p_b(1)))],'r')
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'b');
        patch([nonan_data_plt_b(discon_p_b(i-1)+1:discon_p_b(i)) fliplr(nonan_data_plt_b(discon_p_b(i-1)+1:discon_p_b(i)))], [upper_bound_plt_b(nonan_data_plt_b(discon_p_b(i-1)+1:discon_p_b(i))) lower_bound_plt_b(nonan_data_plt_b(discon_p_b(i-1)+1:discon_p_b(i)))],'r')
        end
        end
        
        title(['GY 2002 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on;
        alpha(0.3) %transparency
        set(gca,'fontsize',13)    
        legend('surf','bot')

cd D:\Àå±â»ıÅÂ\Dynamic\KOEM\gy_2002
%% GY plot

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_no3./14 + obs_std_gy_no3./14;
        lower_bound_plt = obm_gy_no3./14 - obs_std_gy_no3./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_no3),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 55])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);

        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_nh4./14 + obs_std_gy_nh4./14;
        lower_bound_plt = obm_gy_nh4./14 - obs_std_gy_nh4./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_nh4),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 12])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);

        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_chl + obs_std_gy_chl;
        lower_bound_plt = obm_gy_chl - obs_std_gy_chl;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_chl),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);

        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_temp + obs_std_gy_temp;
        lower_bound_plt = obm_gy_temp - obs_std_gy_temp;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_temp),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_salt + obs_std_gy_salt;
        lower_bound_plt = obm_gy_salt - obs_std_gy_salt;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_salt),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency

        title(['GY 2002 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_no3_b./14 + obs_std_gy_no3_b./14;
        lower_bound_plt = obm_gy_no3_b./14 - obs_std_gy_no3_b./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_no3_b),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 55])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
       
                clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_nh4_b./14 + obs_std_gy_nh4_b./14;
        lower_bound_plt = obm_gy_nh4_b./14 - obs_std_gy_nh4_b./14;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_nh4_b),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 12])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        
                clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_chl_b + obs_std_gy_chl_b;
        lower_bound_plt = obm_gy_chl_b - obs_std_gy_chl_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_chl_b),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);

                clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_temp_b + obs_std_gy_temp_b;
        lower_bound_plt = obm_gy_temp_b - obs_std_gy_temp_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_temp_b),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
       
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(obm_gy_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(4) = length(nonan_data_plt); %end point
        upper_bound_plt = obm_gy_salt_b + obs_std_gy_salt_b;
        lower_bound_plt = obm_gy_salt_b - obs_std_gy_salt_b;
        for i = 1:4
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) lower_bound_plt(nonan_data_plt(1:discon_p(1)))],'r');
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))],'r');
        end
        end 
        plot(obm_gy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt_b),upper_bound_plt,'r-','linew',2);
        plot(1:length(obm_gy_salt_b),lower_bound_plt,'r-','linew',2);
        alpha(0.3) %transparency
        
        title(['GY 2002 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot'),'-dpng')    
        



%% GY plot
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_po4,'b','linew',2);
%         plot(obm_gy_po4./14,'r-','linew',2); 
        xlim([1 365]);
%         plot(1:length(obm_gy_po4),obm_gy_po4./14 + obs_std_gy_po4./14,'m-','linew',2);
%         plot(1:length(obm_gy_po4),obm_gy_po4./14 - obs_std_gy_po4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_gy'),'-dpng')


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obm_gy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3),obm_gy_no3./14 + obs_std_gy_no3./14,'m-','linew',2);
        plot(1:length(obm_gy_no3),obm_gy_no3./14 - obs_std_gy_no3./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        ylim([0 55]);
        grid on;
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);
        plot(obm_gy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4),obm_gy_nh4./14 + obs_std_gy_nh4./14,'m-','linew',2);
        plot(1:length(obm_gy_nh4),obm_gy_nh4./14 - obs_std_gy_nh4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        ylim([0 12]);
        grid on;
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);
        plot(obm_gy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
        plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        ylim([0 15])
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);
        plot(obm_gy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp),obm_gy_temp + obs_std_gy_temp,'m-','linew',2);
        plot(1:length(obm_gy_temp),obm_gy_temp - obs_std_gy_temp,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        plot(obm_gy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt),obm_gy_salt + obs_std_gy_salt,'m-','linew',2);
        plot(1:length(obm_gy_salt),obm_gy_salt - obs_std_gy_salt,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_gy'),'-dpng')            
 %% bot
  fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_po4_b,'b','linew',2);
%         plot(obm_gy_po4_b./14,'r-','linew',2);
        xlim([1 365]);
%         plot(1:length(obm_gy_po4_b),obm_gy_po4_b./14 + obs_std_gy_po4_b./14,'m-','linew',2);
%         plot(1:length(obm_gy_po4_b),obm_gy_po4_b./14 - obs_std_gy_po4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_bot_gy'),'-dpng')
 
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        plot(obm_gy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 + obs_std_gy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 - obs_std_gy_no3_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        ylim([0 55]);
        grid on;
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
        plot(obm_gy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 + obs_std_gy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 - obs_std_gy_nh4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        ylim([0 12])
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        plot(obm_gy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_chl_b),obm_gy_chl_b + obs_std_gy_chl_b,'m-','linew',2);
        plot(1:length(obm_gy_chl_b),obm_gy_chl_b - obs_std_gy_chl_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        ylim([0 15])
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);
        plot(obm_gy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_temp_b),obm_gy_temp_b + obs_std_gy_temp_b,'m-','linew',2);
        plot(1:length(obm_gy_temp_b),obm_gy_temp_b - obs_std_gy_temp_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
        plot(obm_gy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_gy_salt_b),obm_gy_salt_b + obs_std_gy_salt_b,'m-','linew',2);
        plot(1:length(obm_gy_salt_b),obm_gy_salt_b - obs_std_gy_salt_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_gy'),'-dpng')    
        
%%south gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_po4,'b','linew',2); xlim([1 365]);
%         plot(obm_sgy_po4./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_po4),obm_sgy_po4./14 + obs_std_sgy_po4./14,'m-','linew',2);
%         plot(1:length(obm_sgy_po4),obm_sgy_po4./14 - obs_std_sgy_po4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_south_sgy'),'-dpng')


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_no3,'b','linew',2);
        plot(obm_sgy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 + obs_std_sgy_no3./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 - obs_std_sgy_no3./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_nh4,'b','linew',2);
        plot(obm_sgy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 + obs_std_sgy_nh4./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 - obs_std_sgy_nh4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_chl,'b','linew',2);
        plot(obm_sgy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl),obm_sgy_chl + obs_std_sgy_chl,'m-','linew',2);
        plot(1:length(obm_sgy_chl),obm_sgy_chl - obs_std_sgy_chl,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_temp,'b','linew',2);
        plot(obm_sgy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp),obm_sgy_temp + obs_std_sgy_temp,'m-','linew',2);
        plot(1:length(obm_sgy_temp),obm_sgy_temp - obs_std_sgy_temp,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_salt,'b','linew',2);
        plot(obm_sgy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt),obm_sgy_salt + obs_std_sgy_salt,'m-','linew',2);
        plot(1:length(obm_sgy_salt),obm_sgy_salt - obs_std_sgy_salt,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_south_sgy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_po4_b,'b','linew',2);xlim([1 365]);
%         plot(obm_sgy_po4_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_sgy_po4_b),obm_sgy_po4_b./14 + obs_std_sgy_po4_b./14,'m-','linew',2);
%         plot(1:length(obm_sgy_po4_b),obm_sgy_po4_b./14 - obs_std_sgy_po4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_bot_south_sgy'),'-dpng')
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_no3_b,'b','linew',2);
        plot(obm_sgy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 + obs_std_sgy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 - obs_std_sgy_no3_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_nh4_b,'b','linew',2);
        plot(obm_sgy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 + obs_std_sgy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 - obs_std_sgy_nh4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_chl_b,'b','linew',2);
        plot(obm_sgy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b + obs_std_sgy_chl_b,'m-','linew',2);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b - obs_std_sgy_chl_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_temp_b,'b','linew',2);
        plot(obm_sgy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b + obs_std_sgy_temp_b,'m-','linew',2);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b - obs_std_sgy_temp_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(sgy_salt_b,'b','linew',2);
        plot(obm_sgy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b + obs_std_sgy_salt_b,'m-','linew',2);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b - obs_std_sgy_salt_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_south_sgy'),'-dpng')    
        
%%east gy plot
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_po4,'b','linew',2); xlim([1 365]);
%         plot(obm_egy_po4./14,'r-','linew',2); 
%         plot(1:length(obm_egy_po4),obm_egy_po4./14 + obs_std_egy_po4./14,'m-','linew',2);
%         plot(1:length(obm_egy_po4),obm_egy_po4./14 - obs_std_egy_po4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_egy'),'-dpng')


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_no3,'b','linew',2);
        plot(obm_egy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 + obs_std_egy_no3./14,'m-','linew',2);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 - obs_std_egy_no3./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_nh4,'b','linew',2);
        plot(obm_egy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 + obs_std_egy_nh4./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 - obs_std_egy_nh4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_chl,'b','linew',2);
        plot(obm_egy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl),obm_egy_chl + obs_std_egy_chl,'m-','linew',2);
        plot(1:length(obm_egy_chl),obm_egy_chl - obs_std_egy_chl,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_temp,'b','linew',2);
        plot(obm_egy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp),obm_egy_temp + obs_std_egy_temp,'m-','linew',2);
        plot(1:length(obm_egy_temp),obm_egy_temp - obs_std_egy_temp,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_salt,'b','linew',2);
        plot(obm_egy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt),obm_egy_salt + obs_std_egy_salt,'m-','linew',2);
        plot(1:length(obm_egy_salt),obm_egy_salt - obs_std_egy_salt,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_egy'),'-dpng')            
 %% bot
 
  fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_po4_b,'b','linew',2);xlim([1 365]);
%         plot(obm_egy_po4_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_egy_po4_b),obm_egy_po4_b./14 + obs_std_egy_po4_b./14,'m-','linew',2);
%         plot(1:length(obm_egy_po4_b),obm_egy_po4_b./14 - obs_std_egy_po4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_bot_egy'),'-dpng')
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_no3_b,'b','linew',2);
        plot(obm_egy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 + obs_std_egy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 - obs_std_egy_no3_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_nh4_b,'b','linew',2);
        plot(obm_egy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 + obs_std_egy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 - obs_std_egy_nh4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_chl_b,'b','linew',2);
        plot(obm_egy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b + obs_std_egy_chl_b,'m-','linew',2);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b - obs_std_egy_chl_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_temp_b,'b','linew',2);
        plot(obm_egy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b + obs_std_egy_temp_b,'m-','linew',2);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b - obs_std_egy_temp_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(egy_salt_b,'b','linew',2);
        plot(obm_egy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b + obs_std_egy_salt_b,'m-','linew',2);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b - obs_std_egy_salt_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_egy'),'-dpng')    
        
%% jinju
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_po4,'b','linew',2);xlim([1 365]);
%         plot(obm_jj_po4./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_po4),obm_jj_po4./14 + obs_std_jj_po4./14,'m-','linew',2);
%         plot(1:length(obm_jj_po4),obm_jj_po4./14 - obs_std_jj_po4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_jj'),'-dpng')



fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_no3,'b','linew',2);
        plot(obm_jj_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 + obs_std_jj_no3./14,'m-','linew',2);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 - obs_std_jj_no3./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_nh4,'b','linew',2);
        plot(obm_jj_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 + obs_std_jj_nh4./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 - obs_std_jj_nh4./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_chl,'b','linew',2);
        plot(obm_jj_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl),obm_jj_chl + obs_std_jj_chl,'m-','linew',2);
        plot(1:length(obm_jj_chl),obm_jj_chl - obs_std_jj_chl,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_temp,'b','linew',2);
        plot(obm_jj_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp),obm_jj_temp + obs_std_jj_temp,'m-','linew',2);
        plot(1:length(obm_jj_temp),obm_jj_temp - obs_std_jj_temp,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_salt,'b','linew',2);
        plot(obm_jj_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt),obm_jj_salt + obs_std_jj_salt,'m-','linew',2);
        plot(1:length(obm_jj_salt),obm_jj_salt - obs_std_jj_salt,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_jj'),'-dpng')            
 %% bot
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_po4_b,'b','linew',2);xlim([1 365]);
%         plot(obm_jj_po4_b./14,'r-','linew',2); xlim([1 365]);
%         plot(1:length(obm_jj_po4_b),obm_jj_po4_b./14 + obs_std_jj_po4_b./14,'m-','linew',2);
%         plot(1:length(obm_jj_po4_b),obm_jj_po4_b./14 - obs_std_jj_po4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL po4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_po4_bot_jj'),'-dpng')
 
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_no3_b,'b','linew',2);
        plot(obm_jj_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 + obs_std_jj_no3_b./14,'m-','linew',2);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 - obs_std_jj_no3_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_no3_bot_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_nh4_b,'b','linew',2);
        plot(obm_jj_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 + obs_std_jj_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 - obs_std_jj_nh4_b./14,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_NH4_bot_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_chl_b,'b','linew',2);
        plot(obm_jj_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b + obs_std_jj_chl_b,'m-','linew',2);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b - obs_std_jj_chl_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_chl_bot_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_temp_b,'b','linew',2);
        plot(obm_jj_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b + obs_std_jj_temp_b,'m-','linew',2);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b - obs_std_jj_temp_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_temp_bot_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(jj_salt_b,'b','linew',2);
        plot(obm_jj_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b + obs_std_jj_salt_b,'m-','linew',2);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b - obs_std_jj_salt_b,'m-','linew',2);
        title(['GY 2002 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('2001_daily_KOEM_OBS_vs_MODEL_salt_bot_jj'),'-dpng')    
        
        
cd D:\Àå±â»ıÅÂ\Dynamic\result\2003
% save('2002_mp_p_sewer_det_result.mat','jj_*', 'gy_*', 'sgy_*', 'egy_*','obm*','obs*');     
save('2003_mp_p_sewer_det_f_result.mat','-v7.3');          
        
        

%plot daily
for i = 1:65
%   for   i = 4:4
   if ind_out(i) == 1
       fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obs_gy_no3(i,:)./14,'r-','linew',2); xlim([1 365]);
        title([name_tag{i},'- 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2002)','fontsize',13)
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
        xlabel('time(days on 2002)','fontsize',13)
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
        xlabel('time(days on 2002)','fontsize',13)
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
        xlabel('time(days on 2002)','fontsize',13)
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
        xlabel('time(days on 2002)','fontsize',13)
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
        xlabel('time(days on 2002)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        name_gf= name_tag{i};
        print(fig,strcat('2001_daily_chl_KOEM_OBS_vs_MODEL_(',num2str(i),')'),'-dpng')
        hold off;
        close;
   end
end
