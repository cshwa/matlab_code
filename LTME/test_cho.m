

close all
clear all
clc

load adcp_20131105.mat
u = SerEmmpersec;
su = AnBTEmmpersec;

[t d] = size(u);

iu = find(u==-32768);
isu = find(su==-32768);

u(iu)=0;
su(isu)=0;
for i = 1:d
    real_u(:,i) = squeeze(u(:,i));% + su;
end

real_u = real_u/10;
% real_u(find(u>=200)) = 0;
real_u(iu) = 0;
contourf(real_u')

h = SerHour;
m = SerMin;
s = SerSec;

total = [h m s real_u];

save u_vel.dat -ascii total