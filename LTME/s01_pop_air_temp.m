% 섬진강 상류 수온의 연간 변화 알아보기
% 환경부 - 구례 자료 사용
% 작성 - 조은별 (2018.08.27)

clc; clear all; close all

%=== 폴더 안에 있는 xls 형식의 파일 연속으로 불러와 원하는 정점만 골라내기===%
folder='D:\mepl\data\03_OtherData\data4lt\database\기상청기온\gy_airtemp\';
t_year = 2015;
%--- 폴더안에 있는 모든 xls 형식의 파일 이름 목록 만들기 ---------------------
f=struct2cell(dir('gy_airtemp/',num2str(t_year),'.xlsx'));

a = size(fic);
m=1;
for i = 1:a  
    filetype=fic(1,i);
    f=fullfile(folder,filetype);
    f = char(f);
    [num,txt,raw] = xlsread(f);
    [b] = size(num,1);
    %--- 기온정보 ------------------------
    % 1. 기온정보
    p_temp = cell2map(raw(2:32,2:13));
  
end

%==========================================================================
