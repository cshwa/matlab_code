clear all; clc; close all;
[dg1,dg2,dg3,tg1,tg2,tg3,pg,tg,cg,sg] = textread('CAST0030.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
[dj1,dj2,dj3,tj1,tj2,tj3,pj,tj,cj,sj] = textread('CAST0103.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);

i = 1 : size(d1,1);
x = find(p(i) > 0.15 & c(i) > 30); 
y = find(p(i) == max(p));

j = [x(1) : y(1)];
k = [y(1)+1 : 89];
down = [p(j) t(j) c(j) s(j)];
up =  [p(k) t(k) c(k) s(k)];

figure;
plot(downg(:,2), downg(:,1),'-k.'); hold on
plot(upg(:,2), upg(:,1),'--k.'); 
plot(downj(:,2), downj(:,1),'-r.');
plot(upj(:,2), upj(:,1),'--r.'); grid on
set(gca,'YDir','reverse');
xlabel('Temp'); ylabel('Pre');

figure;
plot(down(:,4), downg(:,1),'--k.'); hold on
plot(up(:,4), upg(:,1),'-k.'); 
plot(down(:,4), downj(:,1),'--r.');
plot(up(:,4), upj(:,1),'-r.'); grid on
set(gca,'YDir','reverse');
xlabel('Sal'); ylabel('Pre');
