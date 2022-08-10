

close all
clear all
clc

data = load('u_re_time.dat');
t = data(:,1);
dep = data(:,2);
u = data(:,3);

clear data

data = load('s_depth.dat');
sd = data(:,2)-6;
clear data;

fu=[];
for i = 1:3034
    ur = u((i-1)*43+1:i*43);
    deps = dep((i-1)*43+1:i*43);
    ind = find(deps > sd(i));
    ur(ind) = NaN;
    um = nanmean(ur);
    tu = ur - um;
    fu = [fu; tu];
end

data = [t dep fu];
save -ascii u_mean.dat data

uu = zeros(3034,43);
for i = 1:3034
    uu(i,:) = u((i-1)*43+1:i*43);
end
uus = nanmean(uu)';
save -ascii u_mean_all.dat uus