% 1. 같은 날의 위성과 NFRDI자료가 시간대로 있는 파일을 모아 하나의 파일로 만들기
% 2. 이 파일을 읽어 KODC 자료 처리할때처럼 정점별로 resorting
% 3. noaa 와 NFRDI 자료에서 seasonal trend를 제거
% 4. yearly-mean SST trend 비교
% 5. 2,6,8,10 월별 SST trend 비교

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