clear all; clc; close all;
lat = 34 + [57.0 55.2 53.5 51.0 48.3 45.0]./60; lon = 127 + [46.2 48.0 48.0 47.0 48.0 49.0]./60;
i = [1 : 5];
d(i) = spheric_dist(lat(i),lat(i+1),lon(i),lon(i+1));  d = [0 d] ./ 1000; d = cumsum(d);

[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0109.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.15 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+3) : y(1)-2]; d2 = d(2) * ones(size(b,2),1);
st2_2 = sortrows([d2 p(b) t(b) s(b)]);
j = [1 : size(st2_2,1)-1];
k = find(diff(st2_2(j,2)) == 0);
i = 1 : size(k,1);
st2_2(k(i,1),:) = [];

fid = fopen('st2_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st2_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0108.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.15 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+3) : y(1)-2]; d3 = d(3) * ones(size(b,2),1);
st4_2 = sortrows([d3 p(b) t(b) s(b)]);
j = [1 : size(st4_2,1)-1];
k = find(diff(st4_2(j,2)) == 0);
i = 1 : size(k,1);
st4_2(k(i,1),:) = [];

fid = fopen('st4_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st4_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0106.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.15 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+5) : y(1)-20]; d4 = d(4) * ones(size(b,2),1);
st5_2 = sortrows([d4 p(b) t(b) s(b)]); 
j = [1 : size(st5_2,1)-1];
k = find(diff(st5_2(j,2)) == 0);
i = 1 : size(k,1);
st5_2(k(i,1),:) = [];

fid = fopen('st5_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st5_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0105.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.15 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+5) : y(1)-2]; d5 = d(5) * ones(size(b,2),1);
st6_2 = sortrows([d5 p(b) t(b) s(b)]);
j = [1 : size(st6_2,1)-1];
k = find(diff(st6_2(j,2)) == 0);
i = 1 : size(k,1);
st6_2(k(i,1),:) = [];

fid = fopen('st6_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st6_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0104.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.15 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+5) : y(1)]; d6 = d(6) * ones(size(b,2),1);
st7_2 = sortrows([d6 p(b) t(b) s(b)]); 
j = [1 : size(st7_2,1)-1];
k = find(diff(st7_2(j,2)) == 0);
i = 1 : size(k,1);
st7_2(k(i,1),:) = [];

fid = fopen('st7_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st7_2');
fclose(fid);


m12sal = [st2_2;st4_2;st5_2;st6_2;st7_2];
fid = fopen('gysal2.txt','w');
fprintf(fid, '%f %f %f %f\n', m12sal');
fclose(fid);


gsal2 = m12sal;
gsal2(:,2) = -1 .* gsal2(:,2);
fid = fopen('gsal2.txt','w');
fprintf(fid, '%f %f %f %f\n', gsal2');
fclose(fid);


figure; 
plot(st2_2(:,4),st2_2(:,2),'k'); hold on;
plot(st4_2(:,4),st4_2(:,2),'y');
plot(st5_2(:,4),st5_2(:,2),'r');
plot(st6_2(:,4),st6_2(:,2),'g');
plot(st7_2(:,4),st7_2(:,2),'c');
set(gca,'YDir','reverse');
