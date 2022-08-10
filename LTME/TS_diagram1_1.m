% TS?diagram , DO?S and DO?T Diagram
% Input file: E03.dat
% m-file: TS_diagram1.m
% 
clc;clear all;close all
figure();      % T-S diagram �ۼ�
% input_file=input(' Input file  : ','s');   % �̹���� ����ϸ� ���� data�� �ѹ��� �ۼ�
% output_file=input(' Outout file : ','s');
% data_raw=load(input_file);
data_raw=load('E03.dat');    % �ϴ� �� �ڷḸ�� ���� �����ϰ� �Ѵ�.
T=data_raw(:,2);
S=data_raw(:,3);
DO=data_raw(:,4);

P=[0 1000];
dens=20.0:1:30;            % �е� Grid����, ���� �� ����
xy=[32 35 0 27];           % ���� ���� & ���¹���
tsdiagram(S,T,P,dens,'.b',xy);     % tidiagram�� �ۼ��ϱ� ���� ���� (blue ?dot)
print -dpsc TSdiagram1.ps     % TSdiagram1.ps�� ����

% subplot(1,3,2);          % S-DO diagram
% axis('square');          % set plot to "square"
% axis([32 35 3 7]);   % set axis limits
% xlabel('Salinity (psu)');   % place labels on axes
% ylabel('DO (ml/l)');
% hold on;
% plot(S,DO,'Sg','MarkerSize',2);   % green square
% print -dpsc SDOdiagram.ps     % S-DO diagram.ps�� ����

% subplot(1,3,3);      % T-DO diagram
% axis('square');      % set plot to "square"
% axis([0 27 3 7]);    % set axis limits
% ylabel('DO (ml/l)');  % place labels on axes
% xlabel('Temperature ({\circ}C)');
% hold on;
% plot(T,DO,'Dr', 'MarkerSize',2);  % red ������
% print -dpsc TDOdiagram.ps      % T-DO diagram.ps�� ����
%set(gcf,'Color','w');              % ���� White
fclose('all');