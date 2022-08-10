close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
for i = 1997:2015
[raw txt]=xlsread('songjung_discharge_estimate_csh.xlsx',num2str(i),'');

    mth1=raw(1:31, 2);  mth1(isnan(mth1)) = [];           %row=days column=month
    mth2=raw(1:31, 3);  mth2(isnan(mth2)) = [];
    mth3=raw(1:31, 4);  mth3(isnan(mth3)) = [];
    mth4=raw(1:31, 5);  mth4(isnan(mth4)) = [];
    mth5=raw(1:31, 6);  mth5(isnan(mth5)) = [];
    mth6=raw(1:31, 7);  mth6(isnan(mth6)) = [];
    mth7=raw(1:31, 8);  mth7(isnan(mth7)) = [];
    mth8=raw(1:31, 9);  mth8(isnan(mth8)) = [];
    mth9=raw(1:31, 10); mth9(isnan(mth9)) = [];
    mth10=raw(1:31, 11);mth10(isnan(mth10)) = [];
    mth11=raw(1:31, 12);mth11(isnan(mth11)) = [];
    mth12=raw(1:31, 13);mth12(isnan(mth12)) = [];


tt=i-1996;
merg_dis{tt} = [mth1; mth2; mth3; mth4; mth5; mth6; mth7; mth8; mth9; ...
    mth10; mth11; mth12;];
end
pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16}; pre_merg_dis{17}; pre_merg_dis{18}; pre_merg_dis{19};];


interp_t = 1:length(merg_dis);
merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;

save('songjung_estimate_data_1997-2015.mat','pre_merg_dis','merg_dis');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
for i = 1997:2018
[raw txt]=xlsread(strcat('TempExcel송정실시간 일자료',num2str(i),'.xls'),'');
if i ==2018
    mth1=raw(1:31, 2);  mth1(isnan(mth1)) = -999;   %row=days column=month
    mth2=raw(1:31, 3);  mth2(isnan(mth2(1:28))) = -999; mth2(isnan(mth2)) = [];
    mth3=raw(1:31, 4);  mth3(isnan(mth3)) = [];
    mth4=raw(1:31, 5);  mth4(isnan(mth4)) = [];
    mth5=raw(1:31, 6);  mth5(isnan(mth5)) = [];
    mth6=raw(1:31, 7);  mth6(isnan(mth6)) = [];
    mth7=raw(1:31, 8);  mth7(isnan(mth7)) = [];
    mth8=raw(1:31, 9);  mth8(isnan(mth8)) = [];
    mth9=raw(1:31, 10); mth9(isnan(mth9)) = [];
    mth10=raw(1:31, 11);mth10(isnan(mth10)) = [];
    mth11=raw(1:31, 12);mth11(isnan(mth11)) = [];
    mth12=raw(1:31, 13);mth12(isnan(mth12)) = [];
    
else
    mth1=raw(1:31, 2);  mth1(isnan(mth1)) = [];           %row=days column=month
    mth2=raw(1:31, 3);  mth2(isnan(mth2)) = [];
    mth3=raw(1:31, 4);  mth3(isnan(mth3)) = [];
    mth4=raw(1:31, 5);  mth4(isnan(mth4)) = [];
    mth5=raw(1:31, 6);  mth5(isnan(mth5)) = [];
    mth6=raw(1:31, 7);  mth6(isnan(mth6)) = [];
    mth7=raw(1:31, 8);  mth7(isnan(mth7)) = [];
    mth8=raw(1:31, 9);  mth8(isnan(mth8)) = [];
    mth9=raw(1:31, 10); mth9(isnan(mth9)) = [];
    mth10=raw(1:31, 11);mth10(isnan(mth10)) = [];
    mth11=raw(1:31, 12);mth11(isnan(mth11)) = [];
    mth12=raw(1:31, 13);mth12(isnan(mth12)) = [];
end

tt=i-1996;
merg_dis{tt} = [mth1; mth2; mth3; mth4; mth5; mth6; mth7; mth8; mth9; ...
    mth10; mth11; mth12;];
end
pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16}; pre_merg_dis{17}; pre_merg_dis{18}; pre_merg_dis{19}; ...
    pre_merg_dis{20}; pre_merg_dis{21}; pre_merg_dis{22};];

merg_dis(merg_dis == -999) = NaN;

interp_t = 1:length(merg_dis);
merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;

% c_spec = c_spec./255;
% le_yr=2000:2018; % year for legend
% 
% figure; hold on;
% for i = 1:length(merg_dis)
%  plot(merg_dis{i},'color',c_spec(i,:),'linew',1.5);
% end
% hold off;
% legend(num2str(le_yr'));
% alpha(0)

save('songjung_data_1997-2018.mat','pre_merg_dis','merg_dis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\sj_1980to1996
for i = 1980:1996
[raw txt]=xlsread(strcat('songjung_discharge_',num2str(i),'.xls'),'');

mth1=raw(:, 4);  %mth1(isnan(mth1)) = [];          %row=days column=month
%interp NaN
nanx = isnan(mth1);
tdx = 1:numel(mth1); 
mth1(nanx) = interp1(tdx(~nanx), mth1(~nanx), tdx(nanx));

tt=i-1979;
merg_dis{tt} = [mth1;];
clearvars mth1 nanx tdx raw txt
end

% c_spec = [255, 000, 000; 139, 000, 000; 160, 082, 045; 188,143,143; ...
%     255,160,122; 255,165,000; 255,228,181; 255,215,000; 124,252,000; 128,128,000; ...
%     100,000,000; 152,251,152; 000,139,139; 135,206,250; 000,255,255; 000,000,255; ...
%     025,025,112; 211,160,221; 075,000,130; 128,000,139; 255,000,255; 219,112,147; ...
%     192,192,192; 119,136,153; 128,128,128; 047,079,079; 000,000,000; ];
% c_spec = c_spec./255;
c_spec=random_color(length(merg_dis));
le_yr=1980:1996; % year for legend

figure; hold on;
for i = 1:length(merg_dis)
 plot(merg_dis{i},'color',c_spec(i,:),'linew',1.5);
 yr_mean_dis(i) = mean(merg_dis{i});
end
hold off;
legend([num2str(le_yr') repmat('- mean : ',17,1) num2str(yr_mean_dis')]);
alpha(0)

pre_merg_dis=merg_dis;
merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; pre_merg_dis{15}; ...
    pre_merg_dis{16}; pre_merg_dis{17};];
save('songjung_estimate_data_ykang.mat','pre_merg_dis','merg_dis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
for i = 1980:1989
[raw txt]=xlsread('songjung_discharge_estimate_csh.xlsx',num2str(i),'');

mth1=raw(1:31, 2);  mth1(isnan(mth1)) = [];           %row=days column=month
mth2=raw(1:31, 3);  mth2(isnan(mth2)) = [];
mth3=raw(1:31, 4);  mth3(isnan(mth3)) = [];
mth4=raw(1:31, 5);  mth4(isnan(mth4)) = [];
mth5=raw(1:31, 6);  mth5(isnan(mth5)) = [];
mth6=raw(1:31, 7);  mth6(isnan(mth6)) = [];
mth7=raw(1:31, 8);  mth7(isnan(mth7)) = [];
mth8=raw(1:31, 9);  mth8(isnan(mth8)) = [];
mth9=raw(1:31, 10); mth9(isnan(mth9)) = [];
mth10=raw(1:31, 11);mth10(isnan(mth10)) = [];
mth11=raw(1:31, 12);mth11(isnan(mth11)) = [];
mth12=raw(1:31, 13);mth12(isnan(mth12)) = [];

tt=i-1979;
merg_dis{tt} = [mth1; mth2; mth3; mth4; mth5; mth6; mth7; mth8; mth9; ...
    mth10; mth11; mth12;];
end
pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10};];

% interp_t = 1:length(merg_dis);
% merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;
% merg_dis(1)=merg_dis(2); % filling start nan;
save('songjung_estimate_data_csh.mat','pre_merg_dis','merg_dis');

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
load 'songjung_estimate_data_csh.mat'; dis_csh=merg_dis;
load 'songjung_estimate_data_ykang.mat'; dis_ykang=merg_dis;
cd E:\장기생태\Dynamic\06_river\data\댐_수문관련\
load 'gahwa_river_estimate_data_1980s.mat'; dis_gahwa_esti=merg_dis;

%time axe
y = 1980:1989;
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
yr_list{i} = num2str(80+(i-1),'%02d');
end

corr=corrcoef(dis_csh,dis_ykang)
% corr1=corrcoef(dis_dam(1:length(dis_dam_esti)),dis_gahwa);
figure; hold on;  plot(dis_ykang,'linew',1.5); plot(dis_csh,'r','linew',1.5);
title('songjung discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend( 'dis. ykang','dis. csh');
xlabel('Time (year)','fontsize',20);
set(gca,'xtick',t_ax);
set(gca,'xlim',[1 length(dis_csh)]);
set(gca,'xticklabel',yr_list);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
% gtext(strcat('corr (nam & gahwa) = ',num2str(corr1(1,2),2)),'Color','k','FontSize',20)
alpha(0);hold off;


corr=corrcoef(dis_csh,dis_gahwa_esti)
% corr1=corrcoef(dis_dam(1:length(dis_dam_esti)),dis_gahwa);
figure; hold on;  plot(dis_csh,'r','linew',2.5); plot(dis_gahwa_esti,'linew',1.5);
title('gawha & songjung discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('dis. songjung esti', 'dis. gahwa esti');
xlabel('Time (year)','fontsize',20);
set(gca,'xtick',t_ax);
set(gca,'xlim',[1 length(dis_csh)]);
set(gca,'xticklabel',yr_list);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
% gtext(strcat('corr (nam & gahwa) = ',num2str(corr1(1,2),2)),'Color','k','FontSize',20)
alpha(0);hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear ; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
load 'songjung_estimate_data_1997-2015.mat';  dis_song_esti=merg_dis;
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;

%time axe
y = 1997:2015;
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
    if i >= 4
        yr_list{i} = num2str(00+(i-4),'%02d');
    else
        yr_list{i} = num2str(97+(i-1),'%02d');
    end
end

corr=corrcoef(dis_song(1:length(dis_song_esti)),dis_song_esti)
figure; hold on; plot(dis_song(1:length(dis_song_esti)),'linew',2.5); plot(dis_song_esti,'r','linew',1.5); 
title('songjung discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('dis. songjung', 'dis. songjung esti');
xlabel('Time (year)','fontsize',20);
set(gca,'xtick',t_ax);
set(gca,'xlim',[1 length(dis_song_esti)]);
set(gca,'xticklabel',yr_list);
set(gca,'fontsize',13); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
ylim([0 6000])
% gtext(strcat('corr (nam & gahwa) = ',num2str(corr1(1,2),2)),'Color','k','FontSize',20)
alpha(0);hold off;

%%
%%
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
load 'songjung_estimate_data_ykang.mat'; dis_ykang=merg_dis; dis_pre_ykang=pre_merg_dis;
clearvars pre_merg_dis merg_dis
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;

% day on year
y = 1997:2018;
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
for j=1997:2018
for i=1:12
    if i == 1
        pre_t_ax_mth(j-1996,i)=eomday(j,1);
    else
        pre_t_ax_mth(j-1996,i)=sum(eomday(j,1:i));
    end
end
end

for i = 1:length(1997:2018)
    if i == 1
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:);
    else
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:) + sum(t_step(2:i));
    end
end

% make daily mean climate
%%%% Feb.28 == 59 , Feb.29 == 60
sum_dis_present = 0;
for i = 1:length(pre_merg_dis)
temp=pre_merg_dis{i};  %songjung 97~2018

if length(temp) == 366
    temp(60) = [];
end
if i == length(pre_merg_dis)
    temp(temp < 0) = 0; % 2018 has NaN as -999 value
    sum_dis_present = sum_dis_present + temp;
else
    sum_dis_present = sum_dis_present + temp;
end
clearvars temp;
end
clim_song_97_d= sum_dis_present ./ length(pre_merg_dis); %daily climate 1997~2018

plot(clim_song_97_d) 
save('climate_songjung_1997to2018.mat','clim_song_97_d');

%%%%% make  1980s climate
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
load 'songjung_estimate_data_ykang.mat'; dis_ykang=merg_dis; dis_pre_ykang=pre_merg_dis;
clearvars pre_merg_dis merg_dis
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;

% day on year
y = 1980:1989;
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
for j=1980:1989
for i=1:12
    if i == 1
        pre_t_ax_mth(j-1979,i)=eomday(j,1);
    else
        pre_t_ax_mth(j-1979,i)=sum(eomday(j,1:i));
    end
end
end

for i = 1:length(1980:1989)
    if i == 1
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:);
    else
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:) + sum(t_step(2:i));
    end
end

% make daily mean climate
%%%% Feb.28 == 59 , Feb.29 == 60
sum_dis_present = 0;
for i = 1:length(dis_pre_ykang)
temp=dis_pre_ykang{i};  %songjung 97~2018

if length(temp) == 366
    temp(60) = [];
end
    sum_dis_present = sum_dis_present + temp;

clearvars temp;
end
clim_song_80_d= sum_dis_present ./ length(dis_pre_ykang); %daily climate 1997~2018

plot(clim_song_80_d) 
save('climate_songjung_1980to1989.mat','clim_song_80_d');

% plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear; clc; 
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
load climate_songjung_1997to2018.mat
load climate_songjung_1980to1989.mat

figure; plot(clim_song_80_d); hold on;  plot(clim_song_97_d,'r');
mean(clim_song_80_d); sum(clim_song_80_d);
mean(clim_song_97_d); sum(clim_song_97_d);

figure; hold on; plot(clim_song_80_d,'linew',2.5); plot(clim_song_97_d,'r','linew',1.5); 
title('songjung climate discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('80~89', '97~18');
xlabel('Time (month)','fontsize',20);
set(gca,'xtick',pre_t_ax_mth);
set(gca,'xlim',[1 pre_t_ax_mth(end)]);
set(gca,'xticklabel',1:12);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('mean || sum = ',num2str(mean(clim_song_80_d),4),' || ', num2str(sum(clim_song_80_d),3)),'Color','b','FontSize',15)
gtext(strcat('mean || sum  = ',num2str(mean(clim_song_97_d),3),' || ', num2str(sum(clim_song_97_d),3)),'Color','r','FontSize',15)
alpha(0);hold off;


%%%%%%%%%%% make 97to06
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980s
load 'songjung_estimate_data_ykang.mat'; dis_ykang=merg_dis; dis_pre_ykang=pre_merg_dis;
clearvars pre_merg_dis merg_dis
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;

% day on year
y = 1997:2006;
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
for j=1997:2006
for i=1:12
    if i == 1
        pre_t_ax_mth(j-1996,i)=eomday(j,1);
    else
        pre_t_ax_mth(j-1996,i)=sum(eomday(j,1:i));
    end
end
end

for i = 1:length(1997:2006)
    if i == 1
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:);
    else
        pre_t_ax_mth_2(1,1+(i-1)*12:i*12) = pre_t_ax_mth(i,:) + sum(t_step(2:i));
    end
end

% make daily mean climate
%%%% Feb.28 == 59 , Feb.29 == 60
sum_dis_present = 0;
for i = 1:10
temp=pre_merg_dis{i};  %songjung 97~2018

if length(temp) == 366
    temp(60) = [];
end

    sum_dis_present = sum_dis_present + temp;

clearvars temp;
end
clim_song_97to06_d= sum_dis_present ./ 10; %daily climate 1997~2018

plot(clim_song_97to06_d) 
save('climate_songjung_1997to2006.mat','clim_song_97to06_d');

close all; clear; clc; 
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
load climate_songjung_1997to2006.mat
load climate_songjung_1980to1989.mat

figure; hold on; plot(clim_song_80_d,'linew',2.5); plot(clim_song_97to06_d,'r','linew',1.5); 
title('songjung climate discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('80~89', '97~06');
xlabel('Time (month)','fontsize',20);
set(gca,'xtick',pre_t_ax_mth);
set(gca,'xlim',[1 pre_t_ax_mth(end)]);
set(gca,'xticklabel',1:12);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('mean || sum = ',num2str(mean(clim_song_80_d),4),' || ', num2str(sum(clim_song_80_d),3)),'Color','b','FontSize',15)
gtext(strcat('mean || sum  = ',num2str(mean(clim_song_97to06_d),3),' || ', num2str(sum(clim_song_97to06_d),3)),'Color','r','FontSize',15)
alpha(0.1);hold off;



% % make monthly mean
% for i = 1:length(pre_t_ax_mth_2)
% if i == 1
%     dis_mean_mth(i) = mean(dis_song(1:pre_t_ax_mth_2(1)));
% else
%     dis_mean_mth(i) = mean(dis_song(pre_t_ax_mth_2(i-1)+1:pre_t_ax_mth_2(i)));
% end
% end
% for i = 1:12
% clim_song_97_m(i)=mean(dis_mean_mth(i:12:end));%monthly climate 1997~2018
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% make  1980~1991 vs. 1992~2018  juam dem %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980to1996
load 'songjung_estimate_data_ykang.mat'; dis_ykang=merg_dis; dis_pre_ykang=pre_merg_dis;
clearvars pre_merg_dis merg_dis
cd E:\장기생태\Dynamic\06_river\data\sj_00_16_year
load 'songjung_data_1997-2018.mat','pre_merg_dis','merg_dis'; dis_song=merg_dis;

%%% merge discharge
dis_pre_total = dis_pre_ykang;
num_ykang=length(dis_pre_ykang);
for i = 1:length(pre_merg_dis)
    temp = pre_merg_dis{i};
if i == length(pre_merg_dis)
    temp(temp < 0) = NaN; % 2018 has NaN as -999 value
    %interp NaN
    nanx = isnan(temp);
    tdx = 1:numel(temp); 
    temp(nanx) = interp1(tdx(~nanx), temp(~nanx), tdx(nanx));
end
    
dis_pre_total{num_ykang+i} = temp;
clearvars temp;
end

%%% plot
cd E:\장기생태\Dynamic\06_river\data\sj_1980to1996
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

c_spec=random_color(length(dis_pre_total));
le_yr=1980:2018; % year for legend
figure; hold on;
for i = 1:length(dis_pre_total)
 plot(dis_pre_total{i},'color',c_spec(i,:),'linew',1.5);
 yr_mean_dis(i) = mean(dis_pre_total{i});
end
hold off;
legend([num2str(le_yr') repmat('- mean : ',length(dis_pre_total),1) num2str(yr_mean_dis',4)]);
% gtext(strcat('before and after = ',num2str(mean(yr_mean_dis(1:12)),4),' and',blanks(2), num2str(mean(yr_mean_dis(13:end)),4)),'Color','k','FontSize',15)
gtext(['before | after = ',num2str(mean(yr_mean_dis(1:12)),4),' |',blanks(1), num2str(mean(yr_mean_dis(13:end)),4)],'Color','k','FontSize',15)
title('songjung discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
xlabel('Time (month)','fontsize',20);
set(gca,'xtick',pre_t_ax_mth);
set(gca,'xlim',[1 pre_t_ax_mth(end)]);
set(gca,'xticklabel',1:12);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
alpha(0)

save('songjung_discharge_1980to2018.mat', 'dis_pre_total');

%%% climate compare whole 
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980to1996
load songjung_discharge_1980to2018.mat
sum_dis_bf= 0; c=0;
sum_dis_af= 0; c2=0;
for i = 1:length(1980:2018)
temp=dis_pre_total{i};  %songjung 97~2018

if length(temp) == 366
    temp(60) = [];
end
if i <= length(1980:1991) 
    sum_dis_bf = sum_dis_bf + temp;
    c = c+1; %check how many loop
elseif i >= length(1980:1992) 
    sum_dis_af = sum_dis_af + temp;
    c2 = c2+1; %check how many loop
end
clearvars temp;
end
 clim_dis_bf = sum_dis_bf./c;
 clim_dis_af = sum_dis_af./c2;
 
figure; hold on; plot(clim_dis_bf,'linew',2.5); plot(clim_dis_af,'r','linew',1.5); 
title('songjung climate discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('80~91', '92~18');
xlabel('Time (month)','fontsize',20);
set(gca,'xtick',pre_t_ax_mth);
set(gca,'xlim',[1 pre_t_ax_mth(end)]);
set(gca,'xticklabel',1:12);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(['mean = ',num2str(mean(clim_dis_bf),4),' m^3/s'],'Color','b','FontSize',15)
gtext(['mean = ',num2str(mean(clim_dis_af),3),' m^3/s'],'Color','r','FontSize',15)
gtext(['climate diff. = ',num2str(sum(clim_dis_bf.*(3600*24))-sum(clim_dis_af.*(3600*24)),10),' m^3 (ton)'],'Color','k','FontSize',15)
alpha(0.1);hold off;

%%% climate compare for each 12yr 
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\sj_1980to1996
load songjung_discharge_1980to2018.mat
sum_dis_bf= 0; c=0;
sum_dis_af= 0; c2=0;
sum_dis_af2= 0; c3=0;
for i = 1:length(1980:2018)-3
temp=dis_pre_total{i};  %songjung 97~2018

if length(temp) == 366
    temp(60) = [];
end
if i <= length(1980:1991) 
    sum_dis_bf = sum_dis_bf + temp;
    c = c+1; %check how many loop
elseif i >= length(1980:1992) && i <= length(1980:2003)
    sum_dis_af = sum_dis_af + temp;
    c2 = c2+1; %check how many loop
elseif i >= length(1980:2004) 
    sum_dis_af2 = sum_dis_af2 + temp;
    c3 = c3+1; %check how many loop
end
clearvars temp;
end
 clim_dis_bf = sum_dis_bf./c;
 clim_dis_af = sum_dis_af./c2;
 clim_dis_af2 = sum_dis_af2./c3;
 
figure; hold on; plot(clim_dis_bf,'linew',2.5); plot(clim_dis_af,'r','linew',1.5); plot(clim_dis_af2,'g','linew',1.5);
title('songjung climate discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('80~91', '92~03', '04~15');
xlabel('Time (month)','fontsize',20);
set(gca,'xtick',pre_t_ax_mth);
set(gca,'xlim',[1 pre_t_ax_mth(end)]);
set(gca,'xticklabel',1:12);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(['mean = ',num2str(mean(clim_dis_bf),4),' m^3/s'],'Color','b','FontSize',15)
gtext(['mean = ',num2str(mean(clim_dis_af),3),' m^3/s'],'Color','r','FontSize',15)
gtext(['mean = ',num2str(mean(clim_dis_af2),3),' m^3/s'],'Color','g','FontSize',15)
alpha(0.1);hold off;



