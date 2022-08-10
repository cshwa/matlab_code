% m-file: stickdia_exam3.m
% By peter
% ����: x�࿡ Date numberǥ���� �̺���

clc;clear all;close all;
a = textread('file12_Nc.dat','');
Data1 = a(:,6:17);              % ���ϴ� �κи� ���� �� ����ۼ�
Nc = reshape(Data1',[],1)       % ���DATA�� 1�� ���DATA�� ��ȯ

% East Component Extract and make matrix
b = textread('file13_Ec.dat','');
Data2 = b(:,6:17);              % ���ϴ� �κи� ���� �� ����ۼ�
Ec = reshape(Data2',[],1)       % Data2' �� Transpose�� �ǹ�, 1�� �ڷ�� ��ȯ

% date number 
n=length(Ec);  % �ڷ��� ����
dn=1:n;
dn=dn';  % 1�� Data�� ����
FData = [dn,Nc,Ec]                 %Fdata ���ǥ��
fid = fopen('Finaldata.dat','w');
fprintf(fid,'%6d%6d%6d\n',FData');   % FData ���� �״�� ���
fclose(fid)
hold on                 % �� ����ؾ� ��
% plot �Լ��� �̿��� Stick diagram �ۼ�
for i = 1:n;
x1 = dn(i);y1 = 0;
x2 = dn(i)+Ec(i);y2 = Nc(i);
plot([x1 x2],[y1 y2],'b');hold on;
end
hold off
% Label �̸������ϱ� 
title('Stick Diagram using plot function');
xlabel(' Time number');
ylabel('Speed(cm/sec)')
axis equal
axis([0 700 -100 100]);    % AXIS([XMIN XMAX YMIN YMAX])