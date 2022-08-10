close all; close all; clc;

yt = 2018;

if yt == 2019
    cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
    load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;
    cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
    load 'namgang_dam_data_2019.mat'; dis_dam=merg_dis;
    load 'namgang_dam_estimate_data.mat'; dis_dam_esti=merg_dis;
    load sacheon_gate_data_2019.mat; dis_gahwa=merg_dis;
    cd E:\장기생태\Dynamic\06_river\data\댐_수문관련
    load gahwa_river_estimate_data.mat; dis_gahwa_esti=merg_dis;

    %time axe
    y = 2000:2019;
    t = mod(y,4)==0 & (mod(y,100)~=0 | mod(y,400)==0);
    t_step=t.*366; t_step(t==0)=365;
    for i = 1:length(t_step)
        if i ==1 
            pre_t_ax(1)=t_step(1);
        else
            pre_t_ax(i)=sum(t_step(1:i))
        end
    end
    t_ax(1)=1;
    t_ax(2:length(pre_t_ax))=pre_t_ax(1:end-1)+1;
    yr_list={};
    for i = 1:length(t_ax)
    yr_list{i} = num2str(00+(i-1),'%02d');
    end

    corr=corrcoef(dis_dam,dis_gahwa);
    % corr1=corrcoef(dis_dam,dis_gahwa);
    figure; hold on; plot(dis_dam,'r','linew',1.5); plot(dis_gahwa,'linew',1.5); 
    title('namgang-dam vs. gahwa discharge compare (daily mean)','fontsize',20);
    ylabel(gca,'discharge (m^3/s)','fontsize',20);
    legend('dis. Namgang', 'dis. gahwa');
    xlabel('Time (year)','fontsize',20);
    set(gca,'xtick',t_ax);
    set(gca,'xlim',[1 length(dis_dam)]);
    set(gca,'xticklabel',yr_list);
    ylim([0 5000])
    set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
    grid(gca,'on'); 
    gtext(strcat('corr (nam & gahwa) = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
    alpha(0);hold off;

else 
    cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
    load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;
    cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
    load 'namgang_dam_data.mat'; dis_dam=merg_dis;
    load 'namgang_dam_estimate_data.mat'; dis_dam_esti=merg_dis;
    load sacheon_gate_data.mat; dis_gahwa=merg_dis;
    cd E:\장기생태\Dynamic\06_river\data\댐_수문관련
    load gahwa_river_estimate_data.mat; dis_gahwa_esti=merg_dis;

    %time axe
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
    t_ax(1)=1;
    t_ax(2:length(pre_t_ax))=pre_t_ax(1:end-1)+1;
    yr_list={};
    for i = 1:length(t_ax)
    yr_list{i} = num2str(00+(i-1),'%02d');
    end

% % %     plot gahwa vs. songjung
%     corr=corrcoef(dis_song(length(dis_song)-length(dis_gahwa)+1:end),dis_gahwa);
%     % corr1=corrcoef(dis_dam,dis_gahwa);
%     figure; hold on; plot(dis_song(length(dis_song)-length(dis_gahwa)+1:end),'r','linew',1.5); plot(dis_gahwa,'linew',1.5); 
%     title('songjung vs. gahwa discharge compare (daily mean)','fontsize',20);
%     ylabel(gca,'discharge (m^3/s)','fontsize',20);
%     legend('dis. Songjung', 'dis. gahwa');
%     xlabel('Time (year)','fontsize',20);
%     set(gca,'xtick',t_ax);
%     set(gca,'xlim',[1 length(dis_gahwa)]);
%     set(gca,'xticklabel',yr_list);
%     set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
%     grid(gca,'on'); 
%     gtext(strcat('corr (songjung & gahwa) = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
%     alpha(0);hold off;
  

% % %     plot namgang vs. songjung
    corr=corrcoef(dis_song(length(dis_song)-length(dis_gahwa)+1:end),dis_dam);
    % corr1=corrcoef(dis_dam,dis_dam);
    figure; hold on; plot(dis_song(length(dis_song)-length(dis_dam)+1:end),'r','linew',1.5); plot(dis_dam,'linew',1.5); 
    title('songjung vs. Namgang discharge compare (daily mean)','fontsize',20);
    ylabel(gca,'discharge (m^3/s)','fontsize',20);
    legend('dis. Songjung', 'dis. Namgang');
    xlabel('Time (year)','fontsize',20);
    set(gca,'xtick',t_ax);
    set(gca,'xlim',[1 length(dis_dam)]);
    set(gca,'xticklabel',yr_list);
    set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
    grid(gca,'on'); 
    gtext(strcat('corr (songjung & Namgang) = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
    alpha(0);hold off;
    

end

return

%%%%%%%%%%% make 2000 to 2018
close all; close all; clc;
cd D:\장기생태\Dynamic\06_river\data\Namgang_dam
load 'namgang_dam_data.mat','pre_merg_dis','merg_dis';pre_dis_dam=pre_merg_dis; dis_dam=merg_dis;
load 'sacheon_gate_data.mat','pre_merg_dis','merg_dis'; pre_dis_gahwa=pre_merg_dis; dis_gahwa=merg_dis;
clearvars pre_merg_dis merg_dis
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
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
        pre_t_ax_mth(j-1996,i)=eomday(j,1);
    else
        pre_t_ax_mth(j-1996,i)=sum(eomday(j,1:i));
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
    alpha(0);hold off;
    
save('climate_all_2000to2018.mat','clim_*');

%%%%%%%%%%% make 2001 to 2014
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
load 'namgang_dam_data.mat','pre_merg_dis','merg_dis';pre_dis_dam=pre_merg_dis; dis_dam=merg_dis;
load 'sacheon_gate_data.mat','pre_merg_dis','merg_dis'; pre_dis_gahwa=pre_merg_dis; dis_gahwa=merg_dis;
clearvars pre_merg_dis merg_dis
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;

% day on year
y = 2000:2014;
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
for j=2000:2014
for i=1:12
    if i == 1
        pre_t_ax_mth(j-1996,i)=eomday(j,1);
    else
        pre_t_ax_mth(j-1996,i)=sum(eomday(j,1:i));
    end
end
end

for i = 1:length(2000:2014)
    if i == 1
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:);
    else
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:) + sum(t_step(2:i));
    end
end

% make daily mean climate
%%%% Feb.28 == 59 , Feb.29 == 60
sum_dis_present = 0;
for i = 5:length(pre_merg_dis)-4 %songjung 2000~2014
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

clim_song_01to14_d= sum_dis_present ./ length(4:length(pre_merg_dis)); %daily climate 2000~2018

plot(clim_song_01to14_d) 

sum_dis_present = 0;
for i = 2:length(pre_dis_dam)-4 %Namgang 2000~2018
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

clim_nam_01to14_d= sum_dis_present ./ length(4:length(pre_merg_dis)); %daily climate 2000~2018
plot(clim_nam_01to14_d) 

sum_dis_present = 0;
for i = 2:length(pre_dis_gahwa)-4 %Namgang 2000~2018
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

clim_gahwa_01to14_d= sum_dis_present ./ length(4:length(pre_merg_dis)); %daily climate 2000~2018
plot(clim_gahwa_01to14_d) 

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

corr=corrcoef(clim_song_01to14_d,clim_gahwa_01to14_d);
corr1=corrcoef(clim_song_01to14_d,clim_nam_01to14_d);
    figure; hold on; plot(clim_song_01to14_d,'r','linew',1.5); plot(clim_gahwa_01to14_d,'linew',1.5);  plot(clim_nam_01to14_d,'g','linew',1.5);
    title('Climate discharge compare (daily mean)','fontsize',20);
    ylabel(gca,'discharge (m^3/s)','fontsize',20);
    legend('dis. Songjung','dis. gahwa' ,'dis. Namgang');
    xlabel('Time (month)','fontsize',20);
    set(gca,'xtick',pre_t_ax_mth);
    set(gca,'xlim',[1 pre_t_ax_mth(end)]);
    set(gca,'xticklabel',1:12);
    set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
    grid(gca,'on');
    ylim([0 600])
    gtext(strcat('corr (songjung & gahwa) = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
    gtext(strcat('corr (songjung & Namgang) = ',num2str(corr1(1,2),2)),'Color','k','FontSize',20)
    alpha(0);hold off;
    
save('climate_all_2001to2004.mat','clim_*');
