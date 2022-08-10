close all; clear all; clc; 
for i = 2000:2018
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

return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





