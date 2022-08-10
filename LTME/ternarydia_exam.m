% points on it are plotted by this step.
% By Gerry Middleton & peter
% m-file: ternarydia_exam.m
% input data: ternarysample.dat'

% ���� �⺻ �ﰢ���̷��׷� Ʋ�����
clc;clear all;close all
xscale = 2/sqrt(3)
x= [0 xscale/2 xscale 0]; % ��������
y = [-1 0 -1 -1];
plot(x,y,'k','LineWidth',2);   % �ﰢ�� �ܰ�
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

% data�� x,y,z������ ������ 3���� �ڷ�
% �� ������ x,y,z���� ������ 100 �̾�� ��
% input percent���� x�� �ﰢ���̾�׷��� �Ʒ� ���� (SILT), 
% y�� �Ʒ� ������ (CLAY), z�� ���� ������ (SAND)
data = load('ternarysample.dat');
x = data(:,1);y = data(:,2);z = data(:,3);
xscale = 2/sqrt(3);
yy = -(100-z)/100;
x3 = x*xscale/100;
x2 = z/173.2;     % tan 60 is 1.732
xx = xscale - x2 - x3;
plot(xx,yy,'s',...
   'LineWidth',1,...           % ���β�
   'MarkerEdgeColor','b',...  % �Ķ��� �ܰ��� ��
   'MarkerFaceColor','g',...   % ���λ� ���
   'MarkerSize',7)  
text(-0.1, -1,'SILT');text(1.17, -1,'CLAY');text(0.54, 0.04, 'SAND')
set(gcf,'Color','w')   % ���� ������� ����