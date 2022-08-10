clear all; clc; close all;
lat = 34 + [57.0 55.2 53.5 51.0 48.3 45.0]./60; lon = 127 + [46.2 48.0 48.0 47.0 48.0 49.0]./60;
i = [1 : 5];
d(i) = spheric_dist(lat(i),lat(i+1),lon(i),lon(i+1));  d = [0 d] ./ 1000; d = cumsum(d);

[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0096.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [92 : y(1)]; d1 = d(1) * ones(size(b,2),1);
st1_1 = sortrows([d1 p(b) t(b) s(b)]);
j = [1 : size(st1_1,1)-1];
k = find(diff(st1_1(j,2)) == 0);
i = 1 : size(k,1);
st1_1(k(i,1),:) = [];

fid = fopen('st1_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st1_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0097.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.20 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+5) : y(1)]; d2 = d(2) * ones(size(b,2),1);
st2_1 = sortrows([d2 p(b) t(b) s(b)]);
j = [1 : size(st2_1,1)-1];
k = find(diff(st2_1(j,2)) == 0);
i = 1 : size(k,1);
st2_1(k(i,1),:) = [];

fid = fopen('st2_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st2_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0098.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.20 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+5) : y(1)-1]; d3 = d(3) * ones(size(b,2),1);
st4_1 = sortrows([d3 p(b) t(b) s(b)]);
j = [1 : size(st4_1,1)-1];
k = find(diff(st4_1(j,2)) == 0);
i = 1 : size(k,1);
st4_1(k(i,1),:) = []; st4_1(279:281,:) = [];

fid = fopen('st4_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st4_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0100.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.20 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+4) : y(1)-2]; d4 = d(4) * ones(size(b,2),1);
st5_1 = sortrows([d4 p(b) t(b) s(b)]); 
j = [1 : size(st5_1,1)-1];
k = find(diff(st5_1(j,2)) == 0);
i = 1 : size(k,1);
st5_1(k(i,1),:) = [];

fid = fopen('st5_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st5_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0101.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.20 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+2) : y(1)-2]; d5 = d(5) * ones(size(b,2),1);
st6_1 = sortrows([d5 p(b) t(b) s(b)]);
j = [1 : size(st6_1,1)-1];
k = find(diff(st6_1(j,2)) == 0);
i = 1 : size(k,1);
st6_1(k(i,1),:) = [];

fid = fopen('st6_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st6_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0102.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [120 : y(1)-2]; d6 = d(6) * ones(size(b,2),1);
st7_1 = sortrows([d6 p(b) t(b) s(b)]); 
j = [1 : size(st7_1,1)-1];
k = find(diff(st7_1(j,2)) == 0);
i = 1 : size(k,1);
st7_1(k(i,1),:) = [];

fid = fopen('st7_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st7_1');
fclose(fid);


m11sal = [st1_1;st2_1;st4_1;st5_1;st6_1;st7_1];

fid = fopen('gysal1.txt','w');
fprintf(fid, '%f %f %f %f\n', m11sal');
fclose(fid);

gsal1 = m11sal;
gsal1(:,2) = -1 .* gsal1(:,2);
fid = fopen('gsal1.txt','w');
fprintf(fid, '%f %f %f %f\n', gsal1');
fclose(fid);



figure; plot(st1_1(:,4),st1_1(:,2)); hold on;
plot(st2_1(:,4),st2_1(:,2),'k');
plot(st4_1(:,4),st4_1(:,2),'y');
plot(st5_1(:,4),st5_1(:,2),'r');
plot(st6_1(:,4),st6_1(:,2),'g');
plot(st7_1(:,4),st7_1(:,2),'c');
set(gca,'YDir','reverse');
