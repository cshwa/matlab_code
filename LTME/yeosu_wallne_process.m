close all; clear; clc; 
cd D:\장기생태\Dynamic\06_river\하수종말처리장\여수_월내
[raw txt]=xlsread('수질_여천2(여수월내)_일반측정망_199701-200112.xls','수질(일반측정망)','');
[raw2 txt2]=xlsread('수질_여천2(여수월내)_일반측정망_200201-200612.xls','수질(일반측정망)','');
[raw3 txt3]=xlsread('수질_여천2(여수월내)_일반측정망_200701-201112.xls','수질(일반측정망)','');
[raw4 txt4]=xlsread('수질_여천2(여수월내)_일반측정망_201201-201612.xls','수질(일반측정망)','');
% since 2013.07
[raw5 txt5]=xlsread('수질_월내동수로_일반측정망_201301-201712.xls','수질(일반측정망)','');
[raw6 txt6]=xlsread('수질_월내동수로_일반측정망_201801-202112.xls','수질(일반측정망)','');

raw_r =  NaN(size(raw,1),54); raw_r2 =  NaN(size(raw2,1),54); raw_r3 =  NaN(size(raw3,1),54);
raw_r4 =  NaN(size(raw4,1),54); raw_r5 =  NaN(size(raw5,1),54); raw_r6 =  NaN(size(raw6,1),54);

raw_r(1:size(raw,1),1:size(raw,2)) =  raw; raw_r2(1:size(raw2,1),1:size(raw2,2)) = raw2; 
raw_r3(1:size(raw3,1),1:size(raw3,2)) = raw3; raw_r4(1:size(raw4,1),1:size(raw4,2)) = raw4;
raw_r5(1:size(raw5,1),1:size(raw5,2)) = raw5; raw_r6(1:size(raw6,1),1:size(raw6,2)) =raw6;

raw_data = [raw_r; raw_r2; raw_r3; raw_r4; raw_r5; raw_r6;];

temp= raw_data(:,12);
tn= raw_data(:,9).*1000/14;;  %% mg/L to mM/m^3;
tp= raw_data(:,10).*1000./30.973762;  %% mg/L to mM/m^3;
do= raw_data(:,5).*0.7*44.661;  %% mg/L to mM/m^3;
ph= raw_data(:,4);
ss= raw_data(:,8);

% nh4= raw_data(:,29);
% no3= raw_data(:,30);

% colum num = [ ,, , 7, 8, 9,
% {'4 수소이온농도',' 5 용존산소(㎎/L)','6 BOD(㎎/L)','7 COD(㎎/L)','8 부유물질(㎎/L)','9 총질소(T-N)(㎎/L)','10 총인(T-P)(㎎/L)','11 TOC(㎎/L)','12 수온(℃)','13 페놀류(㎎/L)','14 전기전도도(?S/㎝)'
% ,'15 총대장균군수(총대장균군수/100ml)','16 카드뮴(㎎/L)','17 시안(㎎/L)','18 납(㎎/L)','19 6가크롬(㎎/L)','20 비소(㎎/L)','21 수은(㎎/L)','22 구리(㎎/L)','23 아연(㎎/L)','24 크롬(㎎/L)','25 니켈(㎎/L)'
% ,'26 바륨(㎎/L)','27 셀레늄(㎎/L)','28 용존총질소(㎎/L)','29 암모니아성 질소(㎎/L)','30 질산성 질소(㎎/L)','31 용존총인(㎎/L)','32 인산염인(㎎/L)','33 클로로필 a(㎎/㎥)','34 헥사클로로벤젠 (㎍/L)'
% ,'35 분원성대장균군수','36 불소(㎎/L)','37 색도(㎎/L)','38 노말헥산추출물질(㎎/L)','39 용해성망간(㎎/L)','40 용해성철(㎎/L)','음이온계면활성제(㎎/L)','트리클로로에틸렌(㎎/L)','테트라클로로에틸렌(㎎/L)'
% ,'사염화탄소','1.2-디클로로에탄(㎎/L)','디클로로메탄(㎎/L)','벤젠(㎎/L)','폴리크로리네이티트비페닐(㎎/L)','유기인(㎎/L)','안티몬 (㎎/L)','클로로포름 (㎎/L)','디에틸헥실프탈레이트 (㎎/L)','1.4-다이옥세인 (㎎/L)','투명도(m)',
% '유량(㎥/s)','포름알데히드(㎎/L)'}

%% make 1997 to 2018 'yymm' form
k=0
for i = 1997:2020
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

%% make 1997 to 2018 'yymm' form
for j = 1:12
 ref_date_mm{j,1} = [num2str(j,'%02d')];
end

figure; hold on;
plot(tn)
xticklabels(1997:2020);
xticks(1:12:288);
% ylabel(
xlim([0 288]); grid on;
% legend('광양1', '광양2', '광양3', '광양4', '광양5');
xtickangle(45);


figure; hold on;
plot(tp)
xticklabels(1997:2020);
xticks(1:12:288);
xlim([0 288]); grid on;
% legend('광양1', '광양2', '광양3', '광양4', '광양5');
xtickangle(45);

 save('yeosu_wallne_data.mat');