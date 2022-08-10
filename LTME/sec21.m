clear all; clc; close all;
lat = 34 + [55.0 56.2 56.7 57.0 57.0 56.4 55.0 53.5 52.3 51.5]./60; lon = 127 + [50.0 51.2 52.5 53.5 55.5 57.7 58.0 57.5 56.5 55.4]./60;
i = [1 : 9];
d(i) = spheric_dist(lat(i),lat(i+1),lon(i),lon(i+1)); d = [0 d]; d = cumsum(d) ./ 1000;


[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0019.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [115 : y(1)]; d1 = d(1) * ones(size(b,2),1);
st2_1_1 = sortrows([d1 p(b) t(b) s(b)]);
j = [1 : size(st2_1_1,1)-1];
k = find(diff(st2_1_1(j,2)) == 0);
i = 1 : size(k,1);
st2_1_1(k(i,1),:) = [];

fid = fopen('st2_1_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st2_1_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0020.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [73 : y(1)]; d2 = d(2) * ones(size(b,2),1);
st3_1 = sortrows([d2 p(b) t(b) s(b)]);
j = [1 : size(st3_1,1)-1];
k = find(diff(st3_1(j,2)) == 0);
i = 1 : size(k,1);
st3_1(k(i,1),:) = [];

fid = fopen('st3_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st3_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0021.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [42 : y(1)]; d3 = d(3) * ones(size(b,2),1);
st21_1 = sortrows([d3 p(b) t(b) s(b)]);
j = [1 : size(st21_1,1)-1];
k = find(diff(st21_1(j,2)) == 0);
i = 1 : size(k,1);
st21_1(k(i,1),:) = [];

fid = fopen('st21_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st21_1');
fclose(fid);


[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0022.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.15 | c(a) > 40);
y = find(p(a) == max(p));
b = [x(find(diff(x) ~= 1)+7) : y(1)]; d4 = d(4) * ones(size(b,2),1);
st22_1 = sortrows([d4 p(b) t(b) s(b)]); 
j = [1 : size(st22_1,1)-1];
k = find(diff(st22_1(j,2)) == 0);
i = 1 : size(k,1);
st22_1(k(i,1),:) = [];

fid = fopen('st22_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st22_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0023.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);

y = find(p(a) == max(p));
b = [85 : y(1)]; d5 = d(5) * ones(size(b,2),1);
st23_1 = sortrows([d5 p(b) t(b) s(b)]);
j = [1 : size(st23_1,1)-1];
k = find(diff(st23_1(j,2)) == 0);
i = 1 : size(k,1);
st23_1(k(i,1),:) = [];

fid = fopen('st23_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st23_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0024.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [52 : y(1)]; d6 = d(6) * ones(size(b,2),1);
st24_1 = sortrows([d6 p(b) t(b) s(b)]); 
j = [1 : size(st24_1,1)-1];
k = find(diff(st24_1(j,2)) == 0);
i = 1 : size(k,1);
st24_1(k(i,1),:) = [];

fid = fopen('st24_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st24_1');
fclose(fid);


[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0025.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [35 : y(1)]; d6 = d(6) * ones(size(b,2),1);
st2425_1 = sortrows([d6 p(b) t(b) s(b)]); 
j = [1 : size(st2425_1,1)-1];
k = find(diff(st2425_1(j,2)) == 0);
i = 1 : size(k,1);
st2425_1(k(i,1),:) = [];

fid = fopen('st2425_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st2425_1');
fclose(fid);



[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0026.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
x = find(p(a) > 0.20 & c(a) > 40);
y = find(p(a) == max(p));
b = [x(1)+10 : y(1)]; d7 = d(7) * ones(size(b,2),1);
st25_1 = sortrows([d7 p(b) t(b) s(b)]); 
j = [1 : size(st25_1,1)-1];
k = find(diff(st25_1(j,2)) == 0);
i = 1 : size(k,1);
st25_1(k(i,1),:) = [];

fid = fopen('st25_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st25_1');
fclose(fid);


[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0027.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [65 : y(1)]; d8 = d(8) * ones(size(b,2),1);
st26_1 = sortrows([d8 p(b) t(b) s(b)]); 
j = [1 : size(st26_1,1)-1];
k = find(diff(st26_1(j,2)) == 0);
i = 1 : size(k,1);
st26_1(k(i,1),:) = [];

fid = fopen('st26_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st26_1');
fclose(fid);


[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0028.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [75 : y(1)]; d9 = d(9) * ones(size(b,2),1);
st27_1 = sortrows([d9 p(b) t(b) s(b)]); 
j = [1 : size(st27_1,1)-1];
k = find(diff(st27_1(j,2)) == 0);
i = 1 : size(k,1);
st27_1(k(i,1),:) = [];

fid = fopen('st27_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st27_1');
fclose(fid);


[d1,d2,d3,t1,t2,t3,p,t,c,s] = textread('CAST0029.txt','%d-%d-%d %d:%d:%f %f %f %f %f','headerlines',1);
a = 1 : size(d1,1);
y = find(p(a) == max(p));
b = [80 : y(1)]; d10 = d(10) * ones(size(b,2),1);
st28_1 = sortrows([d10 p(b) t(b) s(b)]); 
j = [1 : size(st28_1,1)-1];
k = find(diff(st28_1(j,2)) == 0);
i = 1 : size(k,1);
st28_1(k(i,1),:) = [];

fid = fopen('st28_1.txt','w');
fprintf(fid, '%f %f %f %f\n', st28_1');
fclose(fid);




m21 = [st2_1_1;st3_1;st21_1;st22_1;st23_1;st24_1;st25_1;st26_1;st27_1;st28_1];
fid = fopen('m21.txt','w');
fprintf(fid, '%f %f %f %f\n', m21');
fclose(fid);


m21_1 = [st2_1_1;st3_1;st21_1;st22_1;st23_1;st2425_1;st25_1;st26_1;st27_1;st28_1];
fid = fopen('m21_1.txt','w');
fprintf(fid, '%f %f %f %f\n', m21_1');
fclose(fid);


jtem1 = [st2_1_1(1,1) st3_1(size(st3_1,1),2);st21_1(1,1) st21_1(size(st21_1,1),2);st22_1(1,1) st22_1(size(st22_1,1),2);st23_1(1,1) st23_1(size(st23_1,1),2);st24_1(1,1) st24_1(size(st24_1,1),2);st25_1(1,1) st25_1(size(st25_1,1),2);st26_1(1,1) st26_1(size(st26_1,1),2);st27_1(1,1) st27_1(size(st27_1,1),2);st28_1(1,1) st28_1(size(st28_1,1),2)];
jtem1(:,2) = -1 .* jtem1(:,2);
fid = fopen('jtem1.txt','w');
fprintf(fid, '%f %f\n', jtem1');
fclose(fid);



figure; plot(st2_1_1(:,3),st2_1_1(:,2)); hold on;
plot(st3_1(:,3),st3_1(:,2),'k');
plot(st21_1(:,3),st21_1(:,2),'y');
plot(st22_1(:,3),st22_1(:,2),'r');
plot(st23_1(:,3),st23_1(:,2),'g');
plot(st24_1(:,3),st24_1(:,2),'c');
plot(st25_1(:,3),st25_1(:,2),':b');
plot(st26_1(:,3),st26_1(:,2),'m');
plot(st27_1(:,3),st27_1(:,2),':r');
plot(st28_1(:,3),st28_1(:,2),':k');
set(gca,'YDir','reverse');
