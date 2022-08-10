
close all
clear all
clc

year = 2012;
dis = [];
for i = 1:12
    file = ([num2str(year),num2char(i,2),'.xls']);
%     data = xlsread(file,'c6:c5000');
    [NUMERIC,TXT,RAW]=xlsread(file);
    data = TXT(6:end,3);
    aa = []
    for j = 1:length(data)
        bb = cell2mat(data(j));
        bb = str2num(bb);
        aa = [aa; bb];
    end
    dis = [dis; aa];
end

mean_dis = mean(dis);