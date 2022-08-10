

close all
clear all
clc

[NUMERIC,TXT,RAW]=xlsread('당진-1,2분기조류(2010).xlsx','PC-1(2분기)');

data = NUMERIC;
yy = data(:,2);
mon = data(:,3);
day = data(:,4);
hh = data(:,5);
min = data(:,6);
tem = data(:,7);
dep = data(:,8);
us = data(:,19);
um = data(:,21);
ub = data(:,23);
vs = data(:,25);
vm = data(:,27);
vb = data(:,29);

clear data;

data = [yy mon day hh min dep us vs um vm ub vb tem];
data = data';

ind = (['   yy   mo   dd   hh   mi       dep        us        vs        um        vm        ub        vb       tem']);
fid = fopen('2dangjin.txt','w');
fprintf(fid,'%105s\n',ind);
fprintf(fid,'%5d%5d%5d%5d%5d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',data);
fclose(fid);

ang = -15 * pi/180;
us2 = us*cos(ang) + vs*sin(ang);
vs2 = vs*cos(ang) - us*sin(ang);
um2 = um*cos(ang) + vm*sin(ang);
vm2 = vm*cos(ang) - um*sin(ang);
ub2 = ub*cos(ang) + vb*sin(ang);
vb2 = vb*cos(ang) - ub*sin(ang);

data2 = [yy mon day hh min dep us2 vs2 um2 vm2 ub2 vb2 tem];
data2 = data2';

fid = fopen('2dangjin_rot.txt','w');
fprintf(fid,'%105s\n',ind);
fprintf(fid,'%5d%5d%5d%5d%5d%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n',data2);
fclose(fid);

