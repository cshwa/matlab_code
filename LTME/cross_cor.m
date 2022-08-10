

close all
clear all
clc

data1 = load('uv_sur2.dat');
data2 = load('ele10.txt');

u = data1(20:2180,1);
ele = data2(20:2180,1);
clear data1; clear data2;

[c, lag] = xcorr(u,ele);
c = c';
res = [lag; c];
fid = fopen('res_cross.dat','w');
fprintf(fid,'%10d %20.3f\n',res);
fclose(fid);