% ������ ��� ������ ���� ��ȭ �˾ƺ���
% ȯ��� - ���� �ڷ� ���
% �ۼ� - ������ (2018.08.27)

clc; clear all; close all

%=== ���� �ȿ� �ִ� xls ������ ���� �������� �ҷ��� ���ϴ� ������ ��󳻱�===%
folder='D:\mepl\data\03_OtherData\data4lt\database\���û���\gy_airtemp\';
t_year = 2017;
%--- �����ȿ� �ִ� ��� xls ������ ���� �̸� ��� ����� ---------------------
f = [folder '����' num2str(t_year) '.xlsx'];
[num,txt,raw] = xlsread(f);

ss = num(1:31,1:12);
atemp = reshape(ss,31*12,1);
atemp(isnan(atemp))=[];

% eval(['atemp', num2str(t_year), '.dat atemp -ascii']);
save atemp2017.dat ss -ascii
%==========================================================================
