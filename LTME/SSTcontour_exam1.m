% 이 실행파일은 지도위에 SST 이미지를 삽입하고자 하는 방법임.
% file name=SSTcontour_exam1.m
% colorscale.m이 필요함 (http://www.mathworks.com에서 다운로드)

clc;clear all;close all;clf
m_proj('mercator','lon',[122 135],'lat',[31 39]);%위경도 범위지정
[sn,LAT,LONG,TEMP]=textread('SSTdata.dat','');
X1 = round(min(LONG));X2 = round(max(LONG));
Y1 = round(min(LAT));Y2 = round(max(LAT));
Xp = [X1:0.1:X2];Yp = [Y1:0.1:Y2];
[Xi,Yi] = meshgrid(Xp, Yp);
colormap('jet');
contour_level = 4:2:16;
Zi = griddata(LONG, LAT, TEMP, Xi, Yi);
[XXi YYi]=m_ll2xy(Xi,Yi);
pcolor(XXi,YYi,Zi);hold on;
%shading flat;
shading interp
% 지도관련 옵션
m_gshhs_i('color','k'); % 준고해상도 지도의 해안선 
m_gshhs_i('patch',[1 1 .6]); %지도색을 연노랑으로 설정
m_grid('box','fancy','tickdir','in','linewidth',.5)  %위경도를 표시한다.
m_text(130.5,38.5,'EAST SEA','vertical','top'); %지명을 지도상에 표시한다.
% colorscale.m을 이용한 수직스케일 바 설정
colorscale([0 64],[4 16],2,'vert','position',[0.92 0.25 0.02 0.5]);
set(gca,'YaxisLocation','right');
ylabel('SST (Deg.)')
print -dpsc SST_contour1.ps  % 파일저장
set(gcf,'Color','w')         % 배경색 흰색