clear all; clc; %close all;
gysal2 = textread('gysal2.txt','');
st2 = textread('st2_2.txt',''); pre2 = st2(:,2); salt2 = st2(:,3); salt2 = smooth(pre2,salt2,30,'moving');
st4 = textread('st4_2.txt',''); pre4 = st4(:,2); salt4 = st4(:,3); salt4 = smooth(pre4,salt4,30,'moving');
st5 = textread('st5_2.txt',''); pre5 = st5(:,2); salt5 = st5(:,3); salt5 = smooth(pre5,salt5,30,'moving');
st6 = textread('st6_2.txt',''); pre6 = st6(:,2); salt6 = st6(:,3); salt6 = smooth(pre6,salt6,30,'moving');
st7 = textread('st7_2.txt',''); pre7 = st7(:,2); salt7 = st7(:,3); salt7 = smooth(pre7,salt7,30,'moving');
salt = [salt2; salt4; salt5; salt6; salt7];


dist = gysal2(:,1); pre = gysal2(:,2);
x1 = round(min(dist)); x2 = round(max(dist)); xp = [x1:0.2:x2];
y1 = round(min(pre)); y2 = round(max(pre)); yp = [y1:0.2:y2];
[xi,yi] = meshgrid(xp,yp);
zi = griddata(dist,pre,salt,xi,yi);

figure;
%contourf(xi,yi,zi);
%pcolor(xi,yi,zi);
%surface(xi,yi,zi,'edgecolor','none','facecolor','flat'); 
contour(xi,yi,zi,[31:0.01:33],'fill','on');
set(gca,'FontSize',15,'FontWeight','bold');
hold on;
[C,h] = contour(xi,yi,zi,[31:0.2:33],'Linecolor',[0 0 0]);
clabel(C,h,'FontSize',15,'FontWeight','bold','LabelSpacing',400,'Rotation',0);
%contour(xi,yi,zi,'Linecolor',[1 0 0]);
text(4.311135,-0.5,'st2','FontSize',15,'FontWeight','bold'); text(7.459899,-0.5,'st4','FontSize',15,'FontWeight','bold'); 
text(12.333414,-0.5,'st5','FontSize',15,'FontWeight','bold'); text(17.560412,-0.5,'st6','FontSize',15,'FontWeight','bold'); text(23.859207,-0.5,'st7','FontSize',15,'FontWeight','bold'); 
colorbar('vert'); set(gca,'FontSize',15,'FontWeight','bold');
caxis([31 33]);
ylabel('pre','FontSize',15,'FontWeight','bold'); xlabel('Dist(km)','FontSize',15,'FontWeight','bold');
set(gca,'YDir','reverse'); 