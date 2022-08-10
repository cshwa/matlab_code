% 섬진강 상류 수온의 연간 변화 알아보기
% 환경부 - 구례 자료 사용
% 작성 - 조은별 (2018.08.27)

clc; clear all; close all

%=== 폴더 안에 있는 xls 형식의 파일 연속으로 불러와 원하는 정점만 골라내기===%
folder='D:\mepl\data\03_OtherData\data4lt\database\기상청기온\gy_airtemp\';
t_year = 2017;
%--- 폴더안에 있는 모든 xls 형식의 파일 이름 목록 만들기 ---------------------
f = [folder '광양' num2str(t_year) '.xlsx'];
[num,txt,raw] = xlsread(f);

ss = num(1:31,1:12);
atemp = reshape(ss,31*12,1);
atemp(isnan(atemp))=[];

% eval(['atemp', num2str(t_year), '.dat atemp -ascii']);
save atemp2017.dat ss -ascii
%==========================================================================
