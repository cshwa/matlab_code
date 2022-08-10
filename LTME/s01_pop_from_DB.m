% ���û ��� - ���°� �񱳿�
% �ۼ� - ������ (2018.08.27)
% �ڷ���ó https://data.kma.go.kr/data/grnd/selectAsosRltmList.do?pgmNo=36

clc; clear all; close all

folder='E:\11.���\������_3�ܰ�\2���⵵\data\���û���\gy_airtemp';
fname = ('20180829163624.csv');
f=fullfile(folder,fname);
[num,txt,raw] = xlsread(f);

%-- ���� -------------
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
%-- ���� -------------

%==========================================================================
save namhae_temp.mat ma_temp me_temp mi_temp r_date