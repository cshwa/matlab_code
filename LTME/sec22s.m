clear all; clc; close all;
lat = 34 + [55.0 56.2 56.7 57.0 57.0 56.4 55.0 53.5 52.3 51.5]./60; lon = 127 + [50.0 51.2 52.5 53.5 55.5 57.7 58.0 57.5 56.5 55.4]./60;
i = [1 : 9];
d(i) = spheric_dist(lat(i),lat(i+1),lon(i),lon(i+1)); d = [0 d]; d = cumsum(d) ./ 1000;


[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0042.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [72 : y(1)]; d1 = d(1) * ones(size(b,2),1);
st2_1_2 = sortrows([d1 p(b) t(b) s(b)]);
j = [1 : size(st2_1_2,1)-1];
k = find(diff(st2_1_2(j,2)) == 0);
i = 1 : size(k,1);
st2_1_2(k(i,1),:) = [];

fid = fopen('st2_1_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st2_1_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0041.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [65 : y(1)]; d2 = d(2) * ones(size(b,2),1);
st3_2 = sortrows([d2 p(b) t(b) s(b)]);
j = [1 : size(st3_2,1)-1];
k = find(diff(st3_2(j,2)) == 0);
i = 1 : size(k,1);
st3_2(k(i,1),:) = [];

fid = fopen('st3_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st3_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0038.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [79 : y(1)]; d3 = d(3) * ones(size(b,2),1);
st21_2 = sortrows([d3 p(b) t(b) s(b)]);
j = [1 : size(st21_2,1)-1];
k = find(diff(st21_2(j,2)) == 0);
i = 1 : size(k,1);
st21_2(k(i,1),:) = [];

fid = fopen('st21_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st21_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0037.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [115 : y(1)]; d4 = d(4) * ones(size(b,2),1);
st22_2 = sortrows([d4 p(b) t(b) s(b)]); 
j = [1 : size(st22_2,1)-1];
k = find(diff(st22_2(j,2)) == 0);
i = 1 : size(k,1);
st22_2(k(i,1),:) = [];

fid = fopen('st22_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st22_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0036.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [57 : y(1)]; d5 = d(5) * ones(size(b,2),1);
st23_2 = sortrows([d5 p(b) t(b) s(b)]);
j = [1 : size(st23_2,1)-1];
k = find(diff(st23_2(j,2)) == 0);
i = 1 : size(k,1);
st23_2(k(i,1),:) = [];

fid = fopen('st23_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st23_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0035.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [75 : y(1)]; d6 = d(6) * ones(size(b,2),1);
st24_2 = sortrows([d6 p(b) t(b) s(b)]); 
j = [1 : size(st24_2,1)-1];
k = find(diff(st24_2(j,2)) == 0);
i = 1 : size(k,1);
st24_2(k(i,1),:) = [];

fid = fopen('st24_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st24_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0034.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [74 : y(1)]; d7 = d(7) * ones(size(b,2),1);
st25_2 = sortrows([d7 p(b) t(b) s(b)]); 
j = [1 : size(st25_2,1)-1];
k = find(diff(st25_2(j,2)) == 0);
i = 1 : size(k,1);
st25_2(k(i,1),:) = [];

fid = fopen('st25_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st25_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0033.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [74 : y(1)-1]; d8 = d(8) * ones(size(b,2),1);
st26_2 = sortrows([d8 p(b) t(b) s(b)]); 
j = [1 : size(st26_2,1)-1];
k = find(diff(st26_2(j,2)) == 0);
i = 1 : size(k,1);
st26_2(k(i,1),:) = [];

fid = fopen('st26_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st26_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0032.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [83 : y(1)-15]; d9 = d(9) * ones(size(b,2),1);
st27_2 = sortrows([d9 p(b) t(b) s(b)]); 
j = [1 : size(st27_2,1)-1];
k = find(diff(st27_2(j,2)) == 0);
i = 1 : size(k,1);
st27_2(k(i,1),:) = [];

fid = fopen('st27_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st27_2');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0031.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [92 : y(1)-2]; d10 = d(10) * ones(size(b,2),1);
st28_2 = sortrows([d10 p(b) t(b) s(b)]); 
j = [1 : size(st28_2,1)-1];
k = find(diff(st28_2(j,2)) == 0);
i = 1 : size(k,1);
st28_2(k(i,1),:) = [];

fid = fopen('st28_2.txt','w');
fprintf(fid, '%f %f %f %f\n', st28_2');
fclose(fid);


m22sal2 = [st2_1_2;st3_2;st21_2;st22_2;st23_2;st24_2;st25_2;st26_2;st27_2;st28_2];
fid = fopen('jjsal2.txt','w');
fprintf(fid, '%f %f %f %f\n', m22sal2');
fclose(fid);


jsal2 = m22sal2;
jsal2(:,2) = -1 .* jsal2(:,2);
fid = fopen('jsal2.txt','w');
fprintf(fid, '%f %f %f %f\n', jsal2');
fclose(fid);


figure; plot(st2_1_2(:,4),st2_1_2(:,2)); hold on;
plot(st3_2(:,4),st3_2(:,2),'k');
plot(st21_2(:,4),st21_2(:,2),'y');
plot(st22_2(:,4),st22_2(:,2),'r');
plot(st23_2(:,4),st23_2(:,2),'g');
plot(st24_2(:,4),st24_2(:,2),'c');
plot(st25_2(:,4),st25_2(:,2),':b');
plot(st26_2(:,4),st26_2(:,2),'m');
plot(st27_2(:,4),st27_2(:,2),':r');
plot(st28_2(:,4),st28_2(:,2),':k');
set(gca,'YDir','reverse');
