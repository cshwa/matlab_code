clear all; clc; %close all;
m12 = textread('m12.txt','');
st1 = textread('st1_2.txt',''); pre1 = st1(:,2); temp1 = st1(:,3); temp1 = smooth(pre1,temp1,30,'moving');
st2 = textread('st2_2.txt',''); pre2 = st2(:,2); temp2 = st2(:,3); temp2 = smooth(pre2,temp2,30,'moving');
st4 = textread('st4_2.txt',''); pre4 = st4(:,2); temp4 = st4(:,3); temp4 = smooth(pre4,temp4,30,'moving');
st5 = textread('st5_2.txt',''); pre5 = st5(:,2); temp5 = st5(:,3); temp5 = smooth(pre5,temp5,30,'moving');
st6 = textread('st6_2.txt',''); pre6 = st6(:,2); temp6 = st6(:,3); temp6 = smooth(pre6,temp6,30,'moving');
st7 = textread('st7_2.txt',''); pre7 = st7(:,2); temp7 = st7(:,3); temp7 = smooth(pre7,temp7,30,'moving');
temp = [temp1; temp2; temp4; temp5; temp6; temp7];


dist = m12(:,1); pre = m12(:,2); temp = m12(:,3); sal = m12(:,4);
x1 = round(min(dist)); x2 = round(max(dist)); xp = [x1:0.2:x2];
y1 = round(min(pre)); y2 = round(max(pre)); yp = [y1:0.2:y2];
[xi,yi] = meshgrid(xp,yp);
zi = griddata(dist,pre,temp,xi,yi);

figure;
%contourf(xi,yi,zi);
%pcolor(xi,yi,zi);
%surface(xi,yi,zi,'edgecolor','none','facecolor','flat'); 
contour(xi,yi,zi,[18:0.01:19.5],'fill','on');
set(gca,'FontSize',15,'FontWeight','bold','YDir','reverse');
hold on;
[C,h] = contour(xi,yi,zi,[18:0.2:19],'Linecolor',[0 0 0]);
clabel(C,h,'FontSize',15,'FontWeight','bold','LabelSpacing',400,'Rotation',0);
%contour(xi,yi,zi,'Linecolor',[1 0 0]);
text(0,-0.5,'st1','FontSize',15,'FontWeight','bold'); text(4.311135,-0.5,'st2','FontSize',15,'FontWeight','bold'); text(7.459899,-0.5,'st4','FontSize',15,'FontWeight','bold'); 
text(12.333414,-0.5,'st5','FontSize',15,'FontWeight','bold'); text(17.560412,-0.5,'st6','FontSize',15,'FontWeight','bold'); text(23.859207,-0.5,'st7','FontSize',15,'FontWeight','bold'); 
colorbar('vert'); set(gca,'FontSize',15,'FontWeight','bold');
caxis([17.5 19.5]);
ylabel('pre','FontSize',15,'FontWeight','bold'); xlabel('Dist(km)','FontSize',15,'FontWeight','bold');
set(gca,'YDir','reverse'); 