clear all; clc; close all;
[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0110.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);

i = 1 : size(d1,1);
x = find(p(i) > 0.15 | c(i) > 40); 
y = find(p(i) == max(p));

j = [x(find(diff(x) ~= 1)+1) : y(1)];
k = [y(1)+1 : x(size(x,1))];
down = [p(j) t(j) c(j) s(j)];
up =  [p(k) t(k) c(k) s(k)];

figure;
plot(down(:,2), down(:,1),'-k.'); hold on
plot(up(:,2), up(:,1),'-r.'); grid on
set(gca,'YDir','reverse');
xlabel('Temp'); ylabel('Pre');

figure;
plot(down(:,4), down(:,1),'-k.'); hold on
plot(up(:,4), up(:,1),'-r.'); grid on
set(gca,'YDir','reverse');
xlabel('Sal'); ylabel('Pre');
