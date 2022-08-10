clc;clear all;close all;

filenum=input('how many casts ?','s');
filenum=str2double(filenum);

f_data=[];
for fi=1:1:filenum
    [f, p]=uigetfile('*.*','select cast data file');
        filedir=[p,f];
        disp(filedir);
    data=xlsread(filedir); data(:,1)=fi;

[max_stand_n,max_stand_c]=max(data(:,2));
data=data(1:max_stand_c,:);
data=sortrows(data,1);
f_data=[f_data;data];
end

fid=fopen('timeseris_ctd.dat','w');
fprintf(fid,'%1d %8.4f %8.4f %8.4f\r\n',f_data');
fclose(fid);