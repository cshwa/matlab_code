% ellipse_exam1.m
%
px = 3; % ����
py = 1; % ����
Angle = pi/6; % 30 Degree
% ellipse
deg = [0:0.01:2*pi]
eqx = 'px*sin(deg)';
eqy = 'py*cos(deg)';
xp = eval(eqx);yp = eval(eqy);
plot(xp,yp,'r');hold on
% rotation
rotx = eval('xp.*cos(Angle)-yp.*sin(Angle)');
roty = eval('xp.*sin(Angle)+yp.*cos(Angle)');
plot(rotx,roty,'k');axis equal
xlim([-4 4]);ylim([-4 4]);

title('Plot�Լ��� Ȱ���� Ÿ��ü �ۼ�')
