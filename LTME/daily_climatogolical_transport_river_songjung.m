%%%%%%%%%%% make 2000 to 2018
close all; close all; clc;
cd D:\장기생태\Dynamic\06_river\data\Namgang_dam
load 'namgang_dam_data.mat','pre_merg_dis','merg_dis';pre_dis_dam=pre_merg_dis; dis_dam=merg_dis;
load 'sacheon_gate_data.mat','pre_merg_dis','merg_dis'; pre_dis_gahwa=pre_merg_dis; dis_gahwa=merg_dis;
clearvars pre_merg_dis merg_dis
cd D:\장기생태\Dynamic\06_river\data\sj_00_16_year
load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;

% day on year
y = 2000:2018;
t = mod(y,4)==0 & (mod(y,100)~=0 | mod(y,400)==0);
t_step=t.*366; t_step(t==0)=365;
for i = 1:length(t_step)
    if i ==1 
        pre_t_ax(1)=t_step(1);
    else
        pre_t_ax(i)=sum(t_step(1:i))
    end
end

% day on month  
for j=2000:2018
for i=1:12
    if i == 1
        pre_t_ax_mth(j-1999,i)=eomday(j,1);
    else
        pre_t_ax_mth(j-1999,i)=sum(eomday(j,1:i));
    end
end
end

for i = 1:length(2000:2018)
    if i == 1
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:);
    else
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:) + sum(t_step(2:i));
    end
end

% make daily mean climate
%%%% Feb.28 == 59 , Feb.29 == 60
sum_dis_present = 0;
for i = 4:length(pre_merg_dis) %songjung 2000~2018
temp=pre_merg_dis{i}; 
i
if length(temp) == 366
    temp(60) = [];
end
    temp(temp < 0) = NaN; % 2018 has NaN as -999 value
    nanx = isnan(temp);
    tdx = 1:numel(temp); 
    temp(nanx) = interp1(tdx(~nanx), temp(~nanx), tdx(nanx));
    sum_dis_present = sum_dis_present + temp;

end
clearvars temp;

clim_song_00to18_d= sum_dis_present ./ length(4:length(pre_merg_dis)); %daily climate 2000~2018

plot(clim_song_00to18_d) 

sum_dis_present = 0;
for i = 1:length(pre_dis_dam) %Namgang 2000~2018
temp=pre_dis_dam{i}; 
i
if length(temp) == 366
    temp(60) = [];
end
    nanx = isnan(temp);
    tdx = 1:numel(temp); 
    temp(nanx) = interp1(tdx(~nanx), temp(~nanx), tdx(nanx));
    if find(isnan(temp)) == 1
       temp(1)=temp(2);  % filling start nan;
    end
    sum_dis_present = sum_dis_present + temp;
end
clearvars temp;

clim_nam_00to18_d= sum_dis_present ./ length(4:length(pre_merg_dis)); %daily climate 2000~2018
plot(clim_nam_00to18_d) 

sum_dis_present = 0;
for i = 1:length(pre_dis_gahwa) %Namgang 2000~2018
temp=pre_dis_gahwa{i}; 
i
if length(temp) == 366
    temp(60) = [];
end
    nanx = isnan(temp);
    tdx = 1:numel(temp); 
    temp(nanx) = interp1(tdx(~nanx), temp(~nanx), tdx(nanx));
    if find(isnan(temp)) == 1
       temp(1)=temp(2);  % filling start nan;
    end
    sum_dis_present = sum_dis_present + temp;
end
clearvars temp;

clim_gahwa_00to18_d= sum_dis_present ./ length(4:length(pre_merg_dis)); %daily climate 2000~2018
plot(clim_gahwa_00to18_d) 

%%%%% day axe
clearvars pre_t_ax_mth
y = 1981:1981;
t = mod(y,4)==0 & (mod(y,100)~=0 | mod(y,400)==0);
t_step=t.*366; t_step(t==0)=365;
for i = 1:length(t_step)
    if i ==1 
        pre_t_ax(1)=t_step(1);
    else
        pre_t_ax(i)=sum(t_step(1:i))
    end
end

% day on month  
for j=1981:1981
for i=1:12
    if i == 1
        pre_t_ax_mth(i)=eomday(j,1);
    else
        pre_t_ax_mth(i)=sum(eomday(j,1:i));
    end
end
end
%%%%%%%

corr=corrcoef(clim_song_00to18_d,clim_gahwa_00to18_d);
corr1=corrcoef(clim_song_00to18_d,clim_nam_00to18_d);
    figure; hold on; plot(clim_song_00to18_d,'r','linew',1.5); plot(clim_gahwa_00to18_d,'linew',1.5);  plot(clim_nam_00to18_d,'g','linew',1.5);
    title('Climate discharge compare (daily mean)','fontsize',20);
    ylabel(gca,'discharge (m^3/s)','fontsize',20);
    legend('dis. Songjung','dis. gahwa' ,'dis. Namgang');
    xlabel('Time (month)','fontsize',20);
    set(gca,'xtick',pre_t_ax_mth);
    set(gca,'xlim',[1 pre_t_ax_mth(end)]);
    set(gca,'xticklabel',1:12);
    ylim([0 600])
    set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
    grid(gca,'on'); 
    gtext(strcat('corr (songjung & gahwa) = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
    gtext(strcat('corr (songjung & Namgang) = ',num2str(corr1(1,2),2)),'Color','k','FontSize',20)
    plot(clim_song_00to18_d .* 0.7,'m--','linew',1.5);
        plot(clim_song_00to18_d .* 0.6,'m--','linew',1.5);
    legend('dis. Songjung','dis. gahwa' ,'dis. Namgang','70% songjung');
    alpha(0);
     
for j=2000:2018
for i=1:12
    if i == 1
        tt_mon_ax(j-1999,i)=eomday(j,1);
    else
        tt_mon_ax(j-1999,i)=sum(eomday(j,1:i));
    end
end
end    

tt_ax = [1 pre_t_ax+1];    
    for i = 1:19
    tt_5mon_ax(i) =  tt_ax(i)-1 + tt_mon_ax(i,5)    
    tt_6mon_ax(i) =  tt_ax(i)-1 + tt_mon_ax(i,6)
    tt_8mon_ax(i) =  tt_ax(i)-1 + tt_mon_ax(i,8)
    end

figure; hold on;
for i = 1:19
%     xline(tt_6mon_ax(i),'color','r');
%     xline(tt_5mon_ax(i)+1,'color','r');
    clearvars size_x 
    size_x = length(tt_5mon_ax(i)+1:tt_6mon_ax(i));
    area(tt_5mon_ax(i)+1:tt_6mon_ax(i),repmat(max(dis_gahwa),size_x,1),'FaceColor','r','facealpha',.1);
end
plot(dis_gahwa,'linew',2); xlabel('Time','fontsize',20); set(gca,'xtick',tt_ax);
    set(gca,'xlim',[1 tt_ax(end)]);
    set(gca,'xticklabel',2000:2019);
ylim([1 max(dis_gahwa)])
    title('day obs transp. on the gahwa','fontsize',20);
    ylabel(gca,'discharge (m^3/s)','fontsize',20);
    legend('6 month','raw. gawha');
    xlabel('Time','fontsize',20);

    
figure; hold on;
for i = 1:19
%     xline(tt_6mon_ax(i),'color','r');
%     xline(tt_5mon_ax(i)+1,'color','r');
    clearvars size_x 
    size_x = length(tt_5mon_ax(i)+1:tt_6mon_ax(i));
    area(tt_5mon_ax(i)+1:tt_6mon_ax(i),repmat(max(dis_song),size_x,1),'FaceColor','r','facealpha',.1);
end
plot(dis_song,'linew',2); xlabel('Time','fontsize',20); set(gca,'xtick',tt_ax);
    set(gca,'xlim',[1 tt_ax(end)]);
    set(gca,'xticklabel',2000:2019);
ylim([1 max(dis_song)])
    title('day obs transp. on the songjung','fontsize',20);
    ylabel(gca,'discharge (m^3/s)','fontsize',20);
    legend('6 month','raw. songjung');
    xlabel('Time','fontsize',20);    
    
% alpha(0)
cd D:\장기생태\Dynamic\06_river\환경과학원
% nutrients
load songjung_yymm_raw_data_89to19_3sig.mat  %com_* is raw and 3sig data

yymmdd_txt_c(t_indx(132)+1) % 2000.01.01
  
figure; hold on;
for i = 1:19
%     xline(tt_6mon_ax(i),'color','r');
%     xline(tt_5mon_ax(i)+1,'color','r');
    clearvars size_x 
    size_x = length(tt_5mon_ax(i)+1:tt_6mon_ax(i));
    area(tt_5mon_ax(i)+1:tt_6mon_ax(i),repmat(max(dis_song),size_x,1),'FaceColor','r','facealpha',.1);
end
plot(dis_song,'linew',2); xlabel('Time','fontsize',20); set(gca,'xtick',tt_ax);
    set(gca,'xlim',[1 tt_ax(end)]);
    set(gca,'xticklabel',2000:2019);
ylim([1 max(dis_song)])
    title('day obs transp. on the songjung','fontsize',20);
    ylabel(gca,'discharge (m^3/s)','fontsize',20);
    legend('6 month','raw. songjung');
    xlabel('Time','fontsize',20);
yyaxis right
plot(com_nh4(t_indx(132)+1:end),'r*')



plot(dis_song,com_nh4(t_indx(132)+1:t_indx(end)-1),'r.')


