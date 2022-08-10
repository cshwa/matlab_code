clear all; clc; %close all;
jjsal1 = textread('jjsal1.txt','');
st2_1 = textread('st2_1_1.txt',''); pre2_1 = st2_1(:,2); salt2_1 = st2_1(:,4); salt2_1 = smooth(pre2_1,salt2_1,30,'moving');
st3 = textread('st3_1.txt',''); pre3 = st3(:,2); salt3 = st3(:,4); salt3 = smooth(pre3,salt3,30,'moving');
st21 = textread('st21_1.txt',''); pre21 = st21(:,2); salt21 = st21(:,4); salt21 = smooth(pre21,salt21,30,'moving');
st22 = textread('st22_1.txt',''); pre22 = st22(:,2); salt22 = st22(:,4); salt22 = smooth(pre22,salt22,30,'moving');
st23 = textread('st23_1.txt',''); pre23 = st23(:,2); salt23 = st23(:,4); salt23 = smooth(pre23,salt23,30,'moving');
st24 = textread('st24_1.txt',''); pre24 = st24(:,2); salt24 = st24(:,4); salt24 = smooth(pre24,salt24,30,'moving');
st25 = textread('st25_1.txt',''); pre25 = st25(:,2); salt25 = st25(:,4); salt25 = smooth(pre25,salt25,30,'moving');
st26 = textread('st26_1.txt',''); pre26 = st26(:,2); salt26 = st26(:,4); salt26 = smooth(pre26,salt26,30,'moving');
st27 = textread('st27_1.txt',''); pre27 = st27(:,2); salt27 = st27(:,4); salt27 = smooth(pre27,salt27,30,'moving');
st28 = textread('st28_1.txt',''); pre28 = st28(:,2); salt28 = st28(:,4); salt28 = smooth(pre28,salt28,30,'moving');

salt = [salt2_1; salt3; salt21; salt22; salt23; salt24; salt25; salt26; salt27; salt28];

dist = jjsal1(:,1); pre = jjsal1(:,2);
%temp = smooth(pre,temp,5); %temp = smooth(dist,temp);
x1 = round(min(dist)); x2 = round(max(dist)); xp = [x1:0.2:x2];
y1 = round(min(pre)); y2 = round(max(pre)); yp = [y1:0.2:y2];
[xi,yi] = meshgrid(xp,yp);
zi = griddata(dist,pre,salt,xi,yi);


figure;
%contourf(xi,yi,zi);
%pcolor(xi,yi,zi);
%surface(xi,yi,zi,'edgecolor','none','facecolor','flat'); 
contour(xi,yi,zi,[32:0.01:33],'fill','on');
hold on;
[C,h] = contour(xi,yi,zi,[32:0.1:33],'Linecolor',[0 0 0]);
clabel(C,h,'FontSize',15,'FontWeight','bold','LabelSpacing',400,'Rotation',0);
%contour(xi,yi,zi,'Linecolor',[1 0 0]);
colorbar('vert'); set(gca,'FontSize',15,'FontWeight','bold');

text(0,-0.5,'st2-1','FontSize',15,'FontWeight','bold'); text(2.874207,-0.5,'st3','FontSize',15,'FontWeight','bold'); text(5.054513,-0.5,'st21','FontSize',15,'FontWeight','bold');
text(6.671222,-0.5,'st22','FontSize',15,'FontWeight','bold'); text(9.707565,-0.5,'st23','FontSize',15,'FontWeight','bold'); text(13.227772,-0.5,'st24','FontSize',15,'FontWeight','bold');
text(15.860586,-0.5,'st25','FontSize',15,'FontWeight','bold');text(18.740851,-0.5,'st26','FontSize',15,'FontWeight','bold');text(21.433227,-0.5,'st27','FontSize',15,'FontWeight','bold');text(23.667122,-0.5,'st28','FontSize',15,'FontWeight','bold');

colorbar('vert');
caxis([31 33]);
set(gca,'FontSize',15,'FontWeight','bold','YDir','reverse');
ylabel('pre','FontSize',15,'FontWeight','bold'); xlabel('Dist(km)','FontSize',15,'FontWeight','bold');
