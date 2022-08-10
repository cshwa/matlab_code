% �� ���������� �������� SST �̹����� �����ϰ��� �ϴ� �����.
% file name=SSTcontour_exam1.m
% colorscale.m�� �ʿ��� (http://www.mathworks.com���� �ٿ�ε�)

clc;clear all;close all;clf
m_proj('mercator','lon',[122 135],'lat',[31 39]);%���浵 ��������
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
% �������� �ɼ�
m_gshhs_i('color','k'); % �ذ��ػ� ������ �ؾȼ� 
m_gshhs_i('patch',[1 1 .6]); %�������� ��������� ����
m_grid('box','fancy','tickdir','in','linewidth',.5)  %���浵�� ǥ���Ѵ�.
m_text(130.5,38.5,'EAST SEA','vertical','top'); %������ ������ ǥ���Ѵ�.
% colorscale.m�� �̿��� ���������� �� ����
colorscale([0 64],[4 16],2,'vert','position',[0.92 0.25 0.02 0.5]);
set(gca,'YaxisLocation','right');
ylabel('SST (Deg.)')
print -dpsc SST_contour1.ps  % ��������
set(gcf,'Color','w')         % ���� ���