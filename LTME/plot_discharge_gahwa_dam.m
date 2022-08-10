close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\댐_수문관련\
for i = 1980:1989
[raw txt]=xlsread('gahwa_river_discharge(estimate).xlsx',num2str(i),'');

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

c_spec = [255, 000, 000; 139, 000, 000; 160, 082, 045; 188,143,143; ...
    255,160,122; 255,165,000; 255,228,181; 255,215,000; 124,252,000; 128,128,000; ...
    100,000,000; 152,251,152; 000,139,139; 135,206,250; 000,255,255; 000,000,255; ...
    025,025,112; 211,160,221; 075,000,130; 128,000,139; 255,000,255; 219,112,147; ...
    192,192,192; 119,136,153; 128,128,128; 047,079,079; 000,000,000;];
c_spec = c_spec./255;
le_yr=1980:1989; % year for legend

figure; hold on;
for i = 1:length(merg_dis)
 plot(merg_dis{i},'color',c_spec(i,:),'linew',1.5);
end
hold off;
legend(num2str(le_yr'));
alpha(0)

pre_merg_dis=merg_dis;
merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10};];

save('gahwa_river_estimate_data_1980s.mat','pre_merg_dis','merg_dis');


%% 2000~2015
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\댐_수문관련\
for i = 2000:2015
[raw txt]=xlsread('gahwa_river_discharge(estimate).xlsx',num2str(i),'');

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

tt=i-1999;
merg_dis{tt} = [mth1; mth2; mth3; mth4; mth5; mth6; mth7; mth8; mth9; ...
    mth10; mth11; mth12;];
end

c_spec = [255, 000, 000; 139, 000, 000; 160, 082, 045; 188,143,143; ...
    255,160,122; 255,165,000; 255,228,181; 255,215,000; 124,252,000; 128,128,000; ...
    100,000,000; 152,251,152; 000,139,139; 135,206,250; 000,255,255; 000,000,255; ...
    025,025,112; 211,160,221; 075,000,130; 128,000,139; 255,000,255; 219,112,147; ...
    192,192,192; 119,136,153; 128,128,128; 047,079,079; 000,000,000;];
c_spec = c_spec./255;
le_yr=2000:2018; % year for legend

figure; hold on;
for i = 1:length(merg_dis)
 plot(merg_dis{i},'color',c_spec(i,:),'linew',1.5);
end
hold off;
legend(num2str(le_yr'));
alpha(0)

pre_merg_dis=merg_dis;
merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16};];

save('gahwa_river_estimate_data.mat','pre_merg_dis','merg_dis');
%%%%%%%%%%%%%%% sacheon regulating gate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
for i = 2000:2018
[raw txt]=xlsread('namgang_dam.xlsx',num2str(i),'');

pre_mth1=raw(:, 11); % mth1(isnan(mth1)) = [];           %row=days column=month
mth1=flipud(pre_mth1);

tt=i-1999;
pre_merg_dis{tt} = [mth1];
end
% pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16}; pre_merg_dis{17}; pre_merg_dis{18}; pre_merg_dis{19};];

interp_t = 1:length(merg_dis);
merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;
merg_dis(1)=merg_dis(2); % filling start nan;
save('sacheon_gate_data.mat','pre_merg_dis','merg_dis');

%%%%%%%%%%%%%%% sacheon regulating gate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
for i = 2000:2019
[raw txt]=xlsread('namgang_dam.xlsx',num2str(i),'');

pre_mth1=raw(:, 11); % mth1(isnan(mth1)) = [];           %row=days column=month
mth1=flipud(pre_mth1);

tt=i-1999;
pre_merg_dis{tt} = [mth1];
end
% pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16}; pre_merg_dis{17}; pre_merg_dis{18}; pre_merg_dis{19}; pre_merg_dis{20};];

interp_t = 1:length(merg_dis);
merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;
merg_dis(1)=merg_dis(2); % filling start nan;
save('sacheon_gate_data_2019.mat','pre_merg_dis','merg_dis');

%%%%%%%%%%%%%%%%% Namgang dam data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
for i = 2000:2018
[raw txt]=xlsread('namgang_dam.xlsx',num2str(i),'');

pre_mth1=raw(:, 8); % mth1(isnan(mth1)) = [];           %row=days column=month
mth1=flipud(pre_mth1);

tt=i-1999;
pre_merg_dis{tt} = [mth1];
end
% pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16}; pre_merg_dis{17}; pre_merg_dis{18}; pre_merg_dis{19};];

interp_t = 1:length(merg_dis);
merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;
merg_dis(1)=merg_dis(2); % filling start nan;
save('namgang_dam_data.mat','pre_merg_dis','merg_dis');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
for i = 2000:2019
[raw txt]=xlsread('namgang_dam.xlsx',num2str(i),'');

pre_mth1=raw(:, 8); % mth1(isnan(mth1)) = [];           %row=days column=month
mth1=flipud(pre_mth1);

tt=i-1999;
pre_merg_dis{tt} = [mth1];
end
% pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16}; pre_merg_dis{17}; pre_merg_dis{18}; pre_merg_dis{19}; pre_merg_dis{20};];

interp_t = 1:length(merg_dis);
merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;
merg_dis(1)=merg_dis(2); % filling start nan;
save('namgang_dam_data_2019.mat','pre_merg_dis','merg_dis');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clear all; clc; 
cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
for i = 2000:2015
[raw txt]=xlsread('namgang_dam_estimate.xlsx',num2str(i),'');

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

tt=i-1999;
merg_dis{tt} = [mth1; mth2; mth3; mth4; mth5; mth6; mth7; mth8; mth9; ...
    mth10; mth11; mth12;];
end
pre_merg_dis=merg_dis;

merg_dis=[pre_merg_dis{1}; pre_merg_dis{2}; pre_merg_dis{3}; pre_merg_dis{4};...
    pre_merg_dis{5}; pre_merg_dis{6}; pre_merg_dis{7}; pre_merg_dis{8}; pre_merg_dis{9}; ...
    pre_merg_dis{10}; pre_merg_dis{11}; pre_merg_dis{12}; pre_merg_dis{13}; pre_merg_dis{14}; ...
    pre_merg_dis{15}; pre_merg_dis{16};];

% interp_t = 1:length(merg_dis);
% merg_dis(isnan(merg_dis))=interp1(interp_t(~isnan(merg_dis)),merg_dis(~isnan(merg_dis)), interp_t(isnan(merg_dis))); % filling middle nan;
% merg_dis(1)=merg_dis(2); % filling start nan;
save('namgang_dam_estimate_data.mat','pre_merg_dis','merg_dis');

%%%%%%%%%%%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; close all; clc;
cd E:\장기생태\Dynamic\06_river\data\Namgang_dam
load 'namgang_dam_data.mat'; dis_dam=merg_dis;
load 'namgang_dam_estimate_data.mat'; dis_dam_esti=merg_dis;
load sacheon_gate_data.mat; dis_gahwa=merg_dis;
cd E:\장기생태\Dynamic\06_river\data\댐_수문관련
load gahwa_river_estimate_data.mat; dis_gahwa_esti=merg_dis;

%time axe
y = 2000:2015;
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

corr=corrcoef(dis_dam(1:length(dis_dam_esti)),dis_dam_esti);
corr1=corrcoef(dis_dam(1:length(dis_dam_esti)),dis_gahwa_esti);
figure; hold on; plot(dis_dam_esti,'r','linew',1.5); plot(dis_dam(1:length(dis_dam_esti)),'linew',1.5); plot(dis_gahwa_esti,'g','linew',1.5); 
title('namgang-dam discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('dis. estimate', 'dis. dam h', 'dis. gahwa');
xlabel('Time (year)','fontsize',20);
set(gca,'xtick',t_ax);
set(gca,'xlim',[1 length(dis_dam_esti)]);
set(gca,'xticklabel',yr_list);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
gtext(strcat('corr (nam & gahwa) = ',num2str(corr1(1,2),2)),'Color','k','FontSize',20)
alpha(0);hold off;

corr=corrcoef(dis_gahwa(1:length(dis_gahwa_esti)),dis_gahwa_esti);
% corr1=corrcoef(dis_dam(1:length(dis_dam_esti)),dis_gahwa_esti);
figure; hold on; plot(dis_gahwa(1:length(dis_gahwa_esti)),'b','linew',1.5); plot(dis_gahwa_esti,'r','linew',1.5); 
title('gahwa river discharge compare (daily mean)','fontsize',20);
ylabel(gca,'discharge (m^3/s)','fontsize',20);
legend('dis. gahwa', 'dis. gawha estimate');
xlabel('Time (year)','fontsize',20);
set(gca,'xtick',t_ax);
set(gca,'xlim',[1 length(dis_gahwa_esti)]);
set(gca,'xticklabel',yr_list);
set(gca,'fontsize',18); %set(gca,'XTickLabelRotation',45);
grid(gca,'on'); 
gtext(strcat('corr = ',num2str(corr(1,2),2)),'Color','k','FontSize',20)
% gtext(strcat('corr (nam & gahwa) = ',num2str(corr1(1,2),2)),'Color','k','FontSize',20)
alpha(0);hold off;






