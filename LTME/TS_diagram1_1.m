% TS?diagram , DO?S and DO?T Diagram
% Input file: E03.dat
% m-file: TS_diagram1.m
% 
clc;clear all;close all
figure();      % T-S diagram 작성
% input_file=input(' Input file  : ','s');   % 이방법을 사용하면 여러 data를 한번에 작성
% output_file=input(' Outout file : ','s');
% data_raw=load(input_file);
data_raw=load('E03.dat');    % 일단 한 자료만을 직접 선택하게 한다.
T=data_raw(:,2);
S=data_raw(:,3);
DO=data_raw(:,4);

P=[0 1000];
dens=20.0:1:30;            % 밀도 Grid설정, 범위 및 간격
xy=[32 35 0 27];           % 염분 범위 & 수온범위
tsdiagram(S,T,P,dens,'.b',xy);     % tidiagram을 작성하기 위한 내용 (blue ?dot)
print -dpsc TSdiagram1.ps     % TSdiagram1.ps로 저장

% subplot(1,3,2);          % S-DO diagram
% axis('square');          % set plot to "square"
% axis([32 35 3 7]);   % set axis limits
% xlabel('Salinity (psu)');   % place labels on axes
% ylabel('DO (ml/l)');
% hold on;
% plot(S,DO,'Sg','MarkerSize',2);   % green square
% print -dpsc SDOdiagram.ps     % S-DO diagram.ps로 저장

% subplot(1,3,3);      % T-DO diagram
% axis('square');      % set plot to "square"
% axis([0 27 3 7]);    % set axis limits
% ylabel('DO (ml/l)');  % place labels on axes
% xlabel('Temperature ({\circ}C)');
% hold on;
% plot(T,DO,'Dr', 'MarkerSize',2);  % red 마름모
% print -dpsc TDOdiagram.ps      % T-DO diagram.ps로 저장
%set(gcf,'Color','w');              % 배경색 White
fclose('all');