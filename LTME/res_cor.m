
close all
clear all
clc

data = load('yg_rot_del.txt');
sp = data(:,6);
% sp = -1*sp;
dep = data(:,7);

[R, lag] = xcorr(sp, dep);
R = R';
res = [lag; R];

fid = fopen('res_cor.dat','w');
fprintf(fid,'%15.5f%15.5f\n',res);
fclose(fid);