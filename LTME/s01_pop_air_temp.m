% ������ ��� ������ ���� ��ȭ �˾ƺ���
% ȯ��� - ���� �ڷ� ���
% �ۼ� - ������ (2018.08.27)

clc; clear all; close all

%=== ���� �ȿ� �ִ� xls ������ ���� �������� �ҷ��� ���ϴ� ������ ��󳻱�===%
folder='D:\mepl\data\03_OtherData\data4lt\database\���û���\gy_airtemp\';
t_year = 2015;
%--- �����ȿ� �ִ� ��� xls ������ ���� �̸� ��� ����� ---------------------
f=struct2cell(dir('gy_airtemp/',num2str(t_year),'.xlsx'));

a = size(fic);
m=1;
for i = 1:a  
    filetype=fic(1,i);
    f=fullfile(folder,filetype);
    f = char(f);
    [num,txt,raw] = xlsread(f);
    [b] = size(num,1);
    %--- ������� ------------------------
    % 1. �������
    p_temp = cell2map(raw(2:32,2:13));
  
end

%==========================================================================
