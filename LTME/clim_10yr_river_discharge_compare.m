close all; clear; clc;
cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge


t_start = 1997; t_end = 2000
sj_trans_out_1st = dis_pre_total{t_start-1979}; %dis_pre_total
for i=t_start+1:t_end
    sj_trans_out = dis_pre_total{t_start-1979};
if length(sj_trans_out) == 366
    sj_trans_out(60) = [];
end
    sj_trans_out_1st = sj_trans_out_1st + sj_trans_out;
end

daily_clim_1st=sj_trans_out_1st ./ length(t_start:t_end);


t_start = 2001; t_end = 2010
sj_trans_out_2nd = dis_pre_total{t_start-1979}; %dis_pre_total
for i=t_start+1:t_end
    sj_trans_out = dis_pre_total{t_start-1979};
if length(sj_trans_out) == 366
    sj_trans_out(60) = [];
end
    sj_trans_out_2nd = sj_trans_out_2nd + sj_trans_out;
end

daily_clim_2nd=sj_trans_out_2nd ./ length(t_start:t_end);

figure; hold on;
plot(daily_clim_1st,'r','linew',2); plot(daily_clim_2nd,'b','linew',2);
text(50,3000,['1997~2000 : ',num2str(mean(daily_clim_1st),'%0.2f'),' m^3/s'],'color','r','fontweight','bold');
text(50,2500,['2001~2010 : ',num2str(mean(daily_clim_2nd),'%0.2f'),' m^3/s'],'color','b','fontweight','bold');
xlim([1 365]); alpha(0.3); grid on; ylabel('songjung discharge (m^3/s)'); xlabel('days');
set(gca,'fontsize',13);


close all; clear; clc;
cd D:\장기생태\Dynamic\06_river\data\sj_1980to1996
load songjung_discharge_1980to2018.mat  % pre_merg_dis is discharge


t_start = 1991; t_end = 2000
sj_trans_out_1st = dis_pre_total{t_start-1979}; %dis_pre_total
for i=t_start+1:t_end
    sj_trans_out = dis_pre_total{t_start-1979};
if length(sj_trans_out) == 366
    sj_trans_out(60) = [];
end
    sj_trans_out_1st = sj_trans_out_1st + sj_trans_out;
end

daily_clim_1st=sj_trans_out_1st ./ length(t_start:t_end);


t_start = 2001; t_end = 2010
sj_trans_out_2nd = dis_pre_total{t_start-1979}; %dis_pre_total
for i=t_start+1:t_end
    sj_trans_out = dis_pre_total{t_start-1979};
if length(sj_trans_out) == 366
    sj_trans_out(60) = [];
end
    sj_trans_out_2nd = sj_trans_out_2nd + sj_trans_out;
end

daily_clim_2nd=sj_trans_out_2nd ./ length(t_start:t_end);

figure; hold on;
plot(daily_clim_1st,'b','linew',2); plot(daily_clim_2nd,'r','linew',2);
text(50,3000,['1991~2000 : ',num2str(mean(daily_clim_1st),'%0.2f'),' m^3/s'],'color','b','fontweight','bold');
text(50,2500,['2001~2010 : ',num2str(mean(daily_clim_2nd),'%0.2f'),' m^3/s'],'color','r','fontweight','bold');
xlim([1 365]); alpha(0.3); grid on; ylabel('섬진강 방류량 (m^3/s)'); xlabel('시간(일)');
set(gca,'fontsize',13); ylim([0 3500]);