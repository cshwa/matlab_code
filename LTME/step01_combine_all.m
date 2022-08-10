% 1. ���� ���� ������ NFRDI�ڷᰡ �ð���� �ִ� ������ ��� �ϳ��� ���Ϸ� �����
% 2. �� ������ �о� KODC �ڷ� ó���Ҷ�ó�� �������� resorting
% 3. noaa �� NFRDI �ڷῡ�� seasonal trend�� ����
% 4. yearly-mean SST trend ��
% 5. 2,6,8,10 ���� SST trend ��

clc; close all; clear all;

dir_to_search = './data/';
txtpattern = fullfile(dir_to_search, '*.txt');
dinfo = dir(txtpattern);
a_data = [];
m = length(dinfo);
for K = 1 : m
  filename = fullfile(dir_to_search, dinfo(K).name);  %just the name
  data = load(filename); %load just this file
  a_data = [a_data; data];  
end

save data/a_data.txt a_data -ascii