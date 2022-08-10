


close all
clear all

load uv2.mat
t1 = [SerMon SerDay SerHour SerMin SerSec SerPG4 AnBTEmmpersec SerEmmpersec];
t2 = [SerMon SerDay SerHour SerMin SerSec SerPG4 AnBTNmmpersec SerNmmpersec];


save u.dat t1 -ascii
save v.dat t2 -ascii