% figure â�� ����ǥ�� �ϱ� 
% fun_disp.m
clc;clear all;close all;
subplot(1,2,1)
x = 0:0.02:2;              
noise = 0.01*randn(size(x)); % add noise 
y = 3*x.*exp(-2*x) + noise;  % ����ǥ��            
plot(x,y,'b.');
xlabel('x');ylabel('y')        
title('Plot of y = 3*x*exp(-2*x) + noise');
	
subplot(1,2,2)
plot(x,y,'r.');             
xlabel('x'); ylabel('y');     %  label ǥ��
title('Plot of {\bf\ity} = 3{\bf\itx}\ite^{-2\itx} + {\itnoise}');  %���ڸ�, ����
set(gcf,'Color',[0.8,0.8,0.7])