close all; clear; clc;
[raw txt]=xlsread('yeosu_monthly_1965-2018.xlsx','sheet1','');
load('Yeosu_tide_1980.mat','elev');
[raw2 txt2]=xlsread('Yeosu_UST_TG.xlsx','sheet1','');
[raw1 txt1]=xlsread('Yeosu_UST_TG.xlsx','sheet2','');
load('elev_to_8761hr_v4.mat','ssh','ano_elev');

% raw2(133,1:2)  % that is 1980.01.01
SSH_80 = raw2(133:144,6);
SSH_69 = raw2(:,6);
SSH_yr=raw1(:,6);

% txt(272,2) % that is 1965.01.01
% txt(320,2)
% txt(452,2) % that is 1980.01.01
msl_p_69=raw(319:906,11);
msl_p_80=raw(451:451+11,11);
airT_69 = raw(319:906,3);
airT_80 = raw(451:451+11,3);


format long
corr=corrcoef(msl_p_80(1:12),SSH_80(1:12));
% [r,lags] = xcorr(msl_p_80(1:12),SSH_80(1:12))
% num2str(r,3)

figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,msl_p_80(1:12),1:12,SSH_80(1:12));
title('1980 Monthly Yeosu MSL pressure(KMA) vs. Elev.(UST)','fontsize',20);
ylabel(hax(1),'Pressure (hPa)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','fontsize',20);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',20);
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

% make monthly climate
for i = 1:12
msl_p_climate(i) = nanmean(msl_p_69(i:12:end));
ssh_climate(i) = mean(SSH_69(i:12:end));
airT_climate(i) = mean(airT_69(i:12:end));
end

corr=corrcoef(msl_p_climate,ssh_climate);
[r,lags] = xcorr(msl_p_climate,ssh_climate,'coeff')
num2str(r,3)
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,msl_p_climate,1:12,ssh_climate);
title('Monthly climate Yeosu MSL pressure(KMA) vs. Elev.(UST)','fontsize',20);
ylabel(hax(1),'Pressure (hPa)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','fontsize',20);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',20);
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;


corr=corrcoef(msl_p_69,SSH_69)
figure;hold on;
[hax,hline1,hline2]=plotyy(1:588,msl_p_69,1:588,SSH_69);
title('Monthly Yeosu MSL pressure(KMA) vs. Elev.(UST)','fontsize',20);
ylabel(hax(1),'Pressure (hPa)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','fontsize',20);
set(hline1,'linew',1.5);set(hline2,'linew',1.5);
xlabel('Time (month)','fontsize',20);
set(hax,'xtick',1:60:588);
set(hax,'xlim',[1 588]);
set(hax,'xticklabel',1:60:588);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

corr=corrcoef(msl_p_69(1:end-1),SSH_69(2:end))
figure;hold on;
[hax,hline1,hline2]=plotyy(1:586,msl_p_69(1:end-2),1:586,SSH_69(3:end));
title('Monthly Yeosu MSL pressure(KMA) vs. Elev.(UST) - lag = 1mth ','fontsize',20);
ylabel(hax(1),'Pressure (hPa)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','color','r','fontsize',20);
set(hline1,'linew',1.5);set(hline2,'color','r','linew',1.5);
xlabel('Time (month)','fontsize',20);
set(hax,'xtick',1:60:588);
set(hax,'xlim',[1 588]);
set(hax,'xticklabel',1:60:588);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;


corr=corrcoef(msl_p_climate(1:end-1),ssh_climate(2:end))
[r,lags] = xcorr(msl_p_climate,ssh_climate,'coeff')
num2str(r,3)
figure;hold on;
[hax,hline1,hline2]=plotyy(1:11,msl_p_climate(1:end-1),1:11,ssh_climate(2:end));
title('Monthly climate Yeosu MSL pressure(KMA) vs. Elev.(UST) -lag = 1mth','fontsize',20);
ylabel(hax(1),'Pressure (hPa)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','fontsize',20);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',20);
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compare to airT

corr=corrcoef(airT_80(1:12),SSH_80(1:12))
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,airT_80(1:12),1:12,SSH_80(1:12));
title('1980 Monthly Yeosu airT (KMA) vs. Elev.(UST)','fontsize',20);
ylabel(hax(1),'Temp (^oC)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','fontsize',20);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',20);
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

corr=corrcoef(airT_69,SSH_69)
figure;hold on;
[hax,hline1,hline2]=plotyy(1:588,airT_69,1:588,SSH_69);
title('1980 Monthly airT (KMA) vs. Elev.(UST)','fontsize',20);
ylabel(hax(1),'Temp (^oC)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','fontsize',20);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',20);
set(hax,'xtick',1:60:588);
set(hax,'xlim',[1 588]);
set(hax,'xticklabel',1:60:588);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

%climate
corr=corrcoef(airT_climate,ssh_climate)
[r,lags] = xcorr(airT_climate,ssh_climate,'coeff')
num2str(r,3)
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,airT_climate,1:12,ssh_climate);
title('Monthly climate airT (KMA) vs. Elev.(UST)','fontsize',20);
ylabel(hax(1),'Temp (^oC)','fontsize',20);
ylabel(hax(2),'Elev. (cm)','fontsize',20);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',20);
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

%%%%%%%%% raw diff vs. atm

diff=(ano_elev(1:length(ssh))./100)-(ssh)'; % unit meters

mth_days = [31 29 31 30 31 30 31 31 30 31 30 30]; %model run 365 but 1980 is 366
for i = 1:12
    if i == 1
        days(i) = mth_days(1)
    else
        days(i) = sum(mth_days(1:i))
    end
end

for i = 1:length(diff)/24
    diff_day(i) = mean(diff((i-1)*24+1:i*24));
    (i-1)*24+1
end

for i = 1:length(ano_elev)/24
    ano_elev_day(i) = mean(ano_elev((i-1)*24+1:i*24));
    (i-1)*24+1
end

for i = 1:length(ssh)/24
    ssh_day(i) = mean(ano_elev((i-1)*24+1:i*24));
end


for i = 1:12
    if i ==1
        diff_mth(i)=mean(diff_day(1:days(i)));
    else
        diff_mth(i)=mean(diff_day(days(i-1)+1:days(i)));
        days(i-1)+1
    end
end

corr=corrcoef(airT_80,diff_mth)
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,airT_80,1:12,diff_mth);
title('Monthly climate airT (KMA) vs. diff. Elev.','fontsize',20);
ylabel(hax(1),'Temp (^oC)','fontsize',18);
ylabel(hax(2),'Elev. (m)','fontsize',18);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',18);
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

corr=corrcoef(msl_p_80,diff_mth)
figure;hold on;
[hax,hline1,hline2]=plotyy(1:12,msl_p_80,1:12,diff_mth);
title('Monthly climate Yeosu MSL pressure(KMA) vs. diff. Elev.','fontsize',20);
ylabel(hax(1),'Pressure (hPa)','fontsize',18);
ylabel(hax(2),'Elev. (m)','fontsize',18);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',18);
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

corr=corrcoef(airT_80(1:end-1),diff_mth(2:end)) %lag - 1mth
corr=corrcoef(msl_p_80(1:end-1),diff_mth(2:end)) %lag - 1mth
figure;hold on;
[hax,hline1,hline2]=plotyy(1:11,msl_p_80(1:end-1),1:11,diff_mth(2:end));
title('Monthly climate Yeosu MSL pressure(KMA) vs. diff. Elev.','fontsize',20);
ylabel(hax(1),'Pressure (hPa)','fontsize',18);
ylabel(hax(2),'Elev. (m)','fontsize',18);
set(hline1,'linew',3);set(hline2,'linew',3);
xlabel('Time (month)','fontsize',18);
set(hax,'xtick',1:1:12);
set(hax,'xlim',[1 12]);
set(hax,'xticklabel',1:1:12);
set(hax,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

%% confirm raw diff and TG - TG re
close all; clear; clc;
load TG_1980_G8_nondet.mat
load('elev_to_8761hr_v4.mat','ssh','ano_elev');

diff=(ano_elev(1:length(ssh))./100)-(ssh)'; % unit meters

mth_days = [31 29 31 30 31 30 31 31 30 31 30 30]; %model run 365 but 1980 is 366
for i = 1:12
    if i == 1
        days(i) = mth_days(1)
    else
        days(i) = sum(mth_days(1:i))
    end
end

for i = 1:length(diff)/24
    diff_day(i) = mean(diff((i-1)*24+1:i*24));
    (i-1)*24+1
end

for i = 1:12
    if i ==1
        diff_mth(i)=mean(diff_day(1:days(i)));
    else
        diff_mth(i)=mean(diff_day(days(i-1)+1:days(i)));
        days(i-1)+1
    end
end

raw_minus_tide=((ano_elev(1:length(pout_raw))./100)-(pout_raw./3.281)');

corr=corrcoef(diff,raw_minus_tide(1:length(diff))) 
figure; hold on;
plot(raw_minus_tide(1:length(diff)),'r'); plot(diff); 
title('Compare TG - model vs. TG - TG_r_e','fontsize',20,'fontweight','bold');
xlabel('Time(Month)','fontsize',18,'fontweight','bold');
ylabel('Elevation (m)','fontsize',18,'fontweight','bold');
set(gca,'xtick',1:744:length(ano_elev));
set(gca,'xticklabel',[1:12]); grid('on');
xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',18,'fontweight','bold');
legend('TG - TG_r_e','TG - model','fontsize',13,'fontweight','bold');
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

plot

figure; plot(ano_elev./100,'b');hold on;
plot(pout_raw./3.281,'r'); plot(((ano_elev(1:length(pout_raw))./100)-(pout_raw./3.281)'),'c'); 




%%%%%%%%%%%%%%%%%% confirm raw diff and TG - TG re %%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc;
load TG_1980_G8_nondet.mat
load('elev_to_8761hr_v4.mat','ssh','ano_elev');

diff=(ano_elev(1:length(ssh))./100)-(ssh)'; % unit meters

mth_days = [31 29 31 30 31 30 31 31 30 31 30 30]; %model run 365 but 1980 is 366
for i = 1:12
    if i == 1
        days(i) = mth_days(1)
    else
        days(i) = sum(mth_days(1:i))
    end
end

for i = 1:length(diff)/24
    diff_day(i) = mean(diff((i-1)*24+1:i*24));
    (i-1)*24+1
end

for i = 1:12
    if i ==1
        diff_mth(i)=mean(diff_day(1:days(i)));
    else
        diff_mth(i)=mean(diff_day(days(i-1)+1:days(i)));
        days(i-1)+1
    end
end

raw_minus_tide=((ano_elev(1:length(pout_raw))./100)-(pout_raw./3.281)');

corr=corrcoef(diff,raw_minus_tide(1:length(diff))) 
figure; hold on;
plot(raw_minus_tide(1:length(diff)),'r'); plot(diff); 
title('Compare TG - model vs. TG - TG_r_e','fontsize',20,'fontweight','bold');
xlabel('Time(Month)','fontsize',18,'fontweight','bold');
ylabel('Elevation (m)','fontsize',18,'fontweight','bold');
set(gca,'xtick',1:744:length(ano_elev));
set(gca,'xticklabel',[1:12]); grid('on');
xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',18,'fontweight','bold');
legend('TG - TG_r_e','TG - model','fontsize',13,'fontweight','bold');
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;

plot(diff-raw_minus_tide(1:length(diff)));

corr=corrcoef(diff,raw_minus_tide(1:length(diff))) 
figure; hold on;
plot(raw_minus_tide(1:length(diff)),'r'); plot(diff); plot(diff-raw_minus_tide(1:length(diff)),'c');
title('Compare TG - model vs. TG - TG_r_e','fontsize',20,'fontweight','bold');
xlabel('Time(Month)','fontsize',18,'fontweight','bold');
ylabel('Elevation (m)','fontsize',18,'fontweight','bold');
set(gca,'xtick',1:744:length(ano_elev));
set(gca,'xticklabel',[1:12]); grid('on');
xlim([1 8761]); % ylim([-0.8 0.8]);
set(gca,'fontsize',18,'fontweight','bold');
legend('TG - TG_r_e','TG - model','diff','fontsize',13,'fontweight','bold');
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
hold off;


diff-raw_minus_tide(1:length(diff))


plot

figure; plot(ano_elev./100,'b');hold on;
plot(pout_raw./3.281,'r'); plot(((ano_elev(1:length(pout_raw))./100)-(pout_raw./3.281)'),'c'); 



