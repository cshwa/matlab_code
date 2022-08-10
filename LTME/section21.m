clear all; clc; %close all;
m21 = textread('m21.txt','');
st2_1 = textread('st2_1_1.txt',''); pre2_1 = st2_1(:,2); temp2_1 = st2_1(:,3); temp2_1 = smooth(pre2_1,temp2_1,30,'moving');
st3 = textread('st3_1.txt',''); pre3 = st3(:,2); temp3 = st3(:,3); temp3 = smooth(pre3,temp3,30,'moving');
st21 = textread('st21_1.txt',''); pre21 = st21(:,2); temp21 = st21(:,3); temp21 = smooth(pre21,temp21,30,'moving');
st22 = textread('st22_1.txt',''); pre22 = st22(:,2); temp22 = st22(:,3); temp22 = smooth(pre22,temp22,30,'moving');
st23 = textread('st23_1.txt',''); pre23 = st23(:,2); temp23 = st23(:,3); temp23 = smooth(pre23,temp23,30,'moving');
st24 = textread('st24_1.txt',''); pre24 = st24(:,2); temp24 = st24(:,3); temp24 = smooth(pre24,temp24,30,'moving');
st25 = textread('st25_1.txt',''); pre25 = st25(:,2); temp25 = st25(:,3); temp25 = smooth(pre25,temp25,30,'moving');
st26 = textread('st26_1.txt',''); pre26 = st26(:,2); temp26 = st26(:,3); temp26 = smooth(pre26,temp26,30,'moving');
st27 = textread('st27_1.txt',''); pre27 = st27(:,2); temp27 = st27(:,3); temp27 = smooth(pre27,temp27,30,'moving');
st28 = textread('st28_1.txt',''); pre28 = st28(:,2); temp28 = st28(:,3); temp28 = smooth(pre28,temp28,30,'moving');

temp = [temp2_1; temp3; temp21; temp22; temp23; temp24; temp25; temp26; temp27; temp28];

dist = m21(:,1); pre = m21(:,2); temp = m21(:,3); sal = m21(:,4);
%temp = smooth(pre,temp,5); %temp = smooth(dist,temp);
x1 = round(min(dist)); x2 = round(max(dist)); xp = [x1:0.2:x2];
y1 = round(min(pre)); y2 = round(max(pre)); yp = [y1:0.2:y2];
[xi,yi] = meshgrid(xp,yp);
zi = griddata(dist,pre,temp,xi,yi);


figure;
%contourf(xi,yi,zi);
%pcolor(xi,yi,zi);
%surface(xi,yi,zi,'edgecolor','none','facecolor','flat'); 
contour(xi,yi,zi,[18:0.01:19],'fill','on');
hold on;
[C,h] = contour(xi,yi,zi,[18:0.2:19],'Linecolor',[0 0 0]);
clabel(C,h,'FontSize',15,'FontWeight','bold','LabelSpacing',400,'Rotation',0);
%contour(xi,yi,zi,'Linecolor',[1 0 0]);
colorbar('vert'); set(gca,'FontSize',15,'FontWeight','bold');

text(0,-0.5,'st2-1','FontSize',15,'FontWeight','bold'); text(2.874207,-0.5,'st3','FontSize',15,'FontWeight','bold'); text(5.054513,-0.5,'st21','FontSize',15,'FontWeight','bold');
text(6.671222,-0.5,'st22','FontSize',15,'FontWeight','bold'); text(9.707565,-0.5,'st23','FontSize',15,'FontWeight','bold'); text(13.227772,-0.5,'st24','FontSize',15,'FontWeight','bold');
text(15.860586,-0.5,'st25','FontSize',15,'FontWeight','bold');text(18.740851,-0.5,'st26','FontSize',15,'FontWeight','bold');text(21.433227,-0.5,'st27','FontSize',15,'FontWeight','bold');text(23.667122,-0.5,'st28','FontSize',15,'FontWeight','bold');

colorbar('vert');
caxis([17.5 19.5]);
set(gca,'FontSize',15,'FontWeight','bold','YDir','reverse');
ylabel('pre','FontSize',15,'FontWeight','bold'); xlabel('Dist(km)','FontSize',15,'FontWeight','bold');
