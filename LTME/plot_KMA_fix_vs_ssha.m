%%%%%%%%%%%%%%%%%%%%%%%%%%%% monthly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
close all; clear all; clc;

time_day = datetime(2017,01,01,00,00,00) + days(0:1:364);
time_n1=datenum(time_day(:));
time_str_day=datestr(time_n1,'mm/dd');

close all; clear all; clc;
list1=dir('*monthly_busan*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);
    if i ==1
        busan_msl_mth_pre=temp_mm.msl_mth;
    else
        busan_msl_mth_pre=busan_msl_mth_pre+temp_mm.msl_mth;
    end
       clear temp_mm temp_mm.msl_mth
end
busan_msl_mth = busan_msl_mth_pre/length(list1);
% plot(1:91,busan_msl_yr,'r'); hold on; line(1:91,mean_val_busan);

list1=dir('*monthly_pohang*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);        
    if i ==1
        pohang_msl_mth_pre=temp_mm.msl_mth;
    else
        pohang_msl_mth_pre=pohang_msl_mth_pre+temp_mm.msl_mth;
    end    
       clear temp_mm temp_mm.msl_mth
end
pohang_msl_mth = pohang_msl_mth_pre/length(list1);
% plot(1:68,pohang_msl_yr,'r'); hold on; line(1:68,mean_val_pohang);

list1=dir('*monthly_gwangju*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);        
    if i ==1
        gwangju_msl_mth_pre=temp_mm.msl_mth;
    else
        gwangju_msl_mth_pre=gwangju_msl_mth_pre+temp_mm.msl_mth;
    end     
       clear temp_mm temp_mm.msl_mth
end
gwangju_msl_mth = gwangju_msl_mth_pre/length(list1);
% plot(1:78,gwangju_msl_yr,'r'); hold on; line(1:78,mean_val_gwangju);

list1=dir('*monthly_seoul*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);        
    if i ==1
        seoul_msl_mth_pre=temp_mm.msl_mth;
    else
        seoul_msl_mth_pre=seoul_msl_mth_pre+temp_mm.msl_mth;
    end      
       clear temp_mm temp_mm.msl_mth
end
seoul_msl_mth = seoul_msl_mth_pre/length(list1);

seoul_17=load('msl_monthly_seoul_2017.mat','msl_mth');
pohang_17=load('msl_monthly_pohang_2017.mat','msl_mth');
busan_17=load('msl_monthly_busan_2017.mat','msl_mth');
gwangju_17=load('msl_monthly_gwangju_2017.mat','msl_mth');
seoul_72=load('msl_monthly_seoul_1972.mat','msl_mth');
pohang_72=load('msl_monthly_pohang_1972.mat','msl_mth');
busan_72=load('msl_monthly_busan_1972.mat','msl_mth');
gwangju_72=load('msl_monthly_gwangju_1972.mat','msl_mth');
seoul_66=load('msl_monthly_seoul_1966.mat','msl_mth');
pohang_66=load('msl_monthly_pohang_1966.mat','msl_mth');
busan_66=load('msl_monthly_busan_1966.mat','msl_mth');
gwangju_66=load('msl_monthly_gwangju_1966.mat','msl_mth');
seoul_09=load('msl_monthly_seoul_2009.mat','msl_mth');
pohang_09=load('msl_monthly_pohang_2009.mat','msl_mth');
busan_09=load('msl_monthly_busan_2009.mat','msl_mth');
gwangju_09=load('msl_monthly_gwangju_2009.mat','msl_mth');

load('tg_ust.mat');
load('tg_ust_inch_ul.mat');

ls_list = dir('ASOS_incheon*');name_file=ls_list.name;
[raw1,txt1]=xlsread(name_file,'','');
incheon_msl= raw1(18:end,3); % 32 - mean sea lv. pressure
incheon_time= txt1(19:end,2); % 32 - mean sea lv. pressure

for i = 1:length(incheon_msl)/12
incheon_msl_yr(i)=mean(incheon_msl((i-1)*12+1:12*i));
disp((i-1)*12+1)
end

ls_list = dir('ASOS_ULSAN*');name_file=ls_list.name;
[raw1,txt1]=xlsread(name_file,'','');
ULSAN_msl= raw1(1:end,3); % 32 - mean sea lv. pressure
ULSAN_time= txt1(2:end,2); % 32 - mean sea lv. pressure

for i = 1:length(ULSAN_msl)/12
ULSAN_msl_yr(i)=mean(ULSAN_msl((i-1)*12+1:12*i));
disp((i-1)*12+1)
end

for i = 1 : 12
ULSAN_msl_mth_clim(i)=mean(ULSAN_msl(i:12:end));
incheon_msl_mth_clim(i)=mean(incheon_msl(i:12:end));
ULSAN_mth_clim(i)=mean(ULSAN_mth(i:12:end));
incheon_mth_clim(i)=mean(incheon_mth(i:12:end));
busan_mth_clim(i)=mean(bu_m(i:12:end,6));
pohang_mth_clim(i)=mean(po_m(i:12:end,6));
end

%%% Busan


%%%%%%%%%%%%%%%%%%%
%%% Pohang
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,pohang_17.msl_mth-pohang_msl_mth,1:12,po_m(end-11:end,6)-pohang_mth_clim');
title('Monthly Pohang MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',16); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,pohang_09.msl_mth-pohang_msl_mth,1:12,po_m(end-95:end-95+11,6)-pohang_mth_clim');
title('Monthly Pohang MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:2:92);
set(hax,'xlim',[1 92]);
set(hax,'xticklabel',1926:2:2017);
set(hax,'fontsize',16); set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

%%% Busan
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,busan_17.msl_mth-busan_msl_mth,1:12,bu_m(end-11:end,6)-busan_mth_clim');
title('Monthly Busan MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',16); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

%%% Ulsan 
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,ULSAN_msl(end-11:end)-ULSAN_msl_mth_clim',1:12,ULSAN_mth(end-11:end,6)-ULSAN_mth_clim');
title('Monthly Ulsan MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',16); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

%%% Incheon
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,incheon_msl(end-11:end)-incheon_msl_mth_clim',1:12,incheon_mth(end-11:end,6)-incheon_mth_clim');
title('Monthly Ulsan MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',16); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;



%%%%%%%%%%%%%%%%%%%%%%%%%%%% Yearly %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc;
list1=dir('*monthly_busan*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);        
    busan_msl_yr(i)=mean(temp_mm.msl_mth);     
       clear temp_mm temp_mm.msl_mth
end
mean_val_busan(1:91) = mean(busan_msl_yr);
% plot(1:91,busan_msl_yr,'r'); hold on; line(1:91,mean_val_busan);

list1=dir('*monthly_pohang*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);        
    pohang_msl_yr(i)=mean(temp_mm.msl_mth);     
       clear temp_mm temp_mm.msl_mth
end
mean_val_pohang(1:69) = mean(pohang_msl_yr);
% plot(1:68,pohang_msl_yr,'r'); hold on; line(1:68,mean_val_pohang);
pohang_msl_yr_fix(1)=pohang_msl_yr(1);
pohang_msl_yr_fix(2) = NaN;
pohang_msl_yr_fix(3:69)=pohang_msl_yr(2:end);

list1=dir('*monthly_gwangju*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);        
    gwangju_msl_yr(i)=mean(temp_mm.msl_mth);     
       clear temp_mm temp_mm.msl_mth
end
mean_val_gwangju(1:78) = mean(gwangju_msl_yr);
% plot(1:78,gwangju_msl_yr,'r'); hold on; line(1:78,mean_val_gwangju);

list1=dir('*monthly_seoul*.mat');
for i = 1:length(list1)  %% for 2017 its time index is not shaking  
    temp_mm = load(list1(i).name);        
    seoul_msl_yr(i)=mean(temp_mm.msl_mth);     
       clear temp_mm temp_mm.msl_mth
end
mean_val_seoul(1:92) = mean(seoul_msl_yr);
% plot(1:88,seoul_msl_yr,'r'); hold on; line(1:88,mean_val_seoul);
%%# make it 1926 to 2017 (except 1950~1953)
seoul_msl_yr_fix(1:24)=seoul_msl_yr(1:24);
seoul_msl_yr_fix(25:28) = NaN;
seoul_msl_yr_fix(29:92)=seoul_msl_yr(25:end);

load('tg_ust.mat');
load('tg_ust_inch_ul.mat');

ls_list = dir('ASOS_incheon*');name_file=ls_list.name;
[raw1,txt1]=xlsread(name_file,'','');
incheon_msl= raw1(18:end,3); % 32 - mean sea lv. pressure
incheon_time= txt1(19:end,2); % 32 - mean sea lv. pressure

for i = 1:length(incheon_msl)/12
incheon_msl_yr(i)=mean(incheon_msl((i-1)*12+1:12*i));
disp((i-1)*12+1)
end

ls_list = dir('ASOS_ULSAN*');name_file=ls_list.name;
[raw1,txt1]=xlsread(name_file,'','');
ULSAN_msl= raw1(1:end,3); % 32 - mean sea lv. pressure
ULSAN_time= txt1(2:end,2); % 32 - mean sea lv. pressure

for i = 1:length(ULSAN_msl)/12
ULSAN_msl_yr(i)=mean(ULSAN_msl((i-1)*12+1:12*i));
disp((i-1)*12+1)
end

figure;hold on;
[hax,hline1,hline2]=plotyy(2:92,busan_msl_yr-mean(busan_msl_yr),92-42:92,bu_a(:,6)-mean(bu_a(:,6)));
title('Yearly Busan MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:2:92);
set(hax,'xlim',[1 92]);
set(hax,'xticklabel',1926:2:2017);
set(hax,'fontsize',16); set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

figure;hold on;
[hax,hline1,hline2]=plotyy(92-71:92,ULSAN_msl_yr-mean(ULSAN_msl_yr),92-40:92,ULSAN_yr-mean(ULSAN_yr));
title('Yearly Ulsan MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:2:92);
set(hax,'xlim',[1 92]);
set(hax,'xticklabel',1926:2:2017);
set(hax,'fontsize',16); set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

figure;hold on;
[hax,hline1,hline2]=plotyy(92-65:92,incheon_msl_yr-mean(incheon_msl_yr),92-39:92,incheon_yr-mean(incheon_yr));
title('Yearly Incheon MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:2:92);
set(hax,'xlim',[1 92]);
set(hax,'xticklabel',1926:2:2017);
set(hax,'fontsize',16); set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

figure;hold on;
[hax,hline1,hline2]=plotyy(24:92,pohang_msl_yr_fix-nanmean(pohang_msl_yr_fix),92-40:92,po_a(:,6)-mean(po_a(:,6)));
title('Yearly Pohang MSL pressure(KMA) vs. SSHA(UST)');
ylabel(hax(1),'Pressure (hPa)');
ylabel(hax(2),'SSHA (cm)');
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (date)');
set(hax,'xtick',1:2:92);
set(hax,'xlim',[1 92]);
set(hax,'xticklabel',1926:2:2017);
set(hax,'fontsize',16); set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

figure;hold on;
le1=plot(2:92,busan_msl_yr-nanmean(ULSAN_msl_yr),'g','linew',2); 
le2=plot(92-71:92,ULSAN_msl_yr-mean(ULSAN_msl_yr),'color','b','linew',2); 
le3=plot(24:92,pohang_msl_yr_fix-nanmean(pohang_msl_yr_fix),'r','linew',2);
le4=plot(92-65:92,incheon_msl_yr-mean(incheon_msl_yr),'k','linew',2);
set(gca,'linew',2);
title('Yearly Mean sea lv. pressure - compare KMA');
xlabel('Time (date)');
ylabel('Pressure(hPa)');
% ylim([1014.5 1018]);
set(gca,'xtick',1:2:92);
set(gca,'xlim',[1 92]);
set(gca,'xticklabel',1926:2:2017);
% set(gca,'ytick',1014.5:0.1:1018);
% set(gca,'yticklabel',1014.5:0.1:1018);
set(gca,'fontsize',16);
led=legend([le1 le2 le3 le4],'Busan','Ulsan','Pohang','Incheon','Orientation','horizontal');
set(gca,'XTickLabelRotation',45); 
led.FontSize = 20;
grid(gca,'on');



figure;hold on;
le1=plot(2:92,busan_msl_yr,'g','linew',2); plot(2:92, mean_val_busan,'--g','linew',2);
% le2=plot(92-42:92,seoul_msl_yr_fix,'color','b','linew',2); plot(1:92, mean_val_seoul,'--b','linew',2);
le3=plot(24:92,pohang_msl_yr_fix,'r','linew',2); plot(24:92, mean_val_pohang,'--r','linew',2);
% le4=plot(15:92,gwangju_msl_yr,'k','linew',2); plot(15:92, mean_val_gwangju,'--k','linew',2);
set(gca,'linew',2);
title('Yearly Mean sea lv. pressure - compare KMA');
xlabel('Time (date)');
ylabel('Pressure(hPa)');
ylim([1014.5 1018]);
set(gca,'xtick',1:2:92);
set(gca,'xticklabel',1926:2:2017);
set(gca,'ytick',1014.5:0.1:1018);
set(gca,'yticklabel',1014.5:0.1:1018);
set(gca,'fontsize',16);
led=legend([le1 le2 le3 le4],'Busan','Seoul','Pohang','Gwangju','Orientation','horizontal');
set(gca,'XTickLabelRotation',45); 
led.FontSize = 20;
grid(gca,'on');


