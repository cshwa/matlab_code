% points on it are plotted by this step.
% By Gerry Middleton & peter
% m-file: ternarydia_exam.m
% input data: ternarysample.dat'

% 먼저 기본 삼각다이러그램 틀만들기
clc;clear all;close all
xscale = 2/sqrt(3)
x= [0 xscale/2 xscale 0]; % 정점결정
y = [-1 0 -1 -1];
plot(x,y,'k','LineWidth',2);   % 삼각형 외각
axis('equal');
axis off;      % turn off rectangular axes
hold on
for i =1:9;    % plot dotted lines parallel to sides,
    x1 = i*xscale/20;  % at 10% interval
    x2 = (20-i)*xscale/20;
    x3 = i*xscale/10;
    y1 = -1 + i/10;y2 = -i/10;
    xx = [x1 x2];yy = [y1 y1];
    xx2 = [x1 x3];yy2 = [y1 -1];
    xx3 = [x3 (xscale/2 + x1)];yy3 = [-1 y2];
 plot (xx,yy,':b',xx2,yy2,':b',xx3,yy3,':b')
end

% data는 x,y,z값으로 구성된 3열의 자료
% 단 각행의 x,y,z값의 총합이 100 이어야 함
% input percent값은 x는 삼각다이어그램의 아래 왼쪽 (SILT), 
% y는 아래 오른쪽 (CLAY), z는 위쪽 꼭지점 (SAND)
data = load('ternarysample.dat');
x = data(:,1);y = data(:,2);z = data(:,3);
xscale = 2/sqrt(3);
yy = -(100-z)/100;
x3 = x*xscale/100;
x2 = z/173.2;     % tan 60 is 1.732
xx = xscale - x2 - x3;
plot(xx,yy,'s',...
   'LineWidth',1,...           % 선두께
   'MarkerEdgeColor','b',...  % 파란색 외곽선 색
   'MarkerFaceColor','g',...   % 연두색 면색
   'MarkerSize',7)  
text(-0.1, -1,'SILT');text(1.17, -1,'CLAY');text(0.54, 0.04, 'SAND')
set(gcf,'Color','w')   % 배경색 흰색으로 지정