% 기상청 기온 - 수온과 비교용
% 작성 - 조은별 (2018.08.27)
% 자료출처 https://data.kma.go.kr/data/grnd/selectAsosRltmList.do?pgmNo=36

clc; clear all; close all

folder='E:\11.사업\장기생태_3단계\2차년도\data\기상청기온\gy_airtemp';
fname = ('20180829163624.csv');
f=fullfile(folder,fname);
[num,txt,raw] = xlsread(f);

%-- 광양 -------------
[m n]= find(num(:,1)>265)
a = min(m);
a=1;
me_temp = num(a:end,3);
ma_temp = num(a:end,6);
mi_temp = num(a:end,4);

date = txt(a+1:end,2);
DateStrings = date;
t = datetime(DateStrings,'InputFormat','yyyy-MM-dd');
r_date = day(t,'dayofyear');

figure()
hold on;
plot(r_date,me_temp,'k.');
plot(r_date,ma_temp,'r.');
plot(r_date,mi_temp,'b.');
%-- 여수 -------------

%==========================================================================
save namhae_temp.mat ma_temp me_temp mi_temp r_date