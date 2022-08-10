% line.m을 이용한 수온 및 염분 수직분포도
% m-file: ts_line.m
% Data load
data = load('CTD_TED6.dat');
depth = data(:,2);
temp = data(:,3);
salinity = data(:,7);
oxygen = data(:,9);

% 수심에 따른 수온분포 그리기
subplot(1,2,1)
LinHd(1)=line(temp, -1*depth);  
set(LinHd, 'color', 'b');
grid on
% label 만들기
xlabel('Temp.(degree)');
ylabel('Depth (m)');
title('TED6 수온 수직분포도');
subplot(1,2,2)
% 수심에 따른 염분분포도
LinHd= line(salinity, -1*depth);  
set(LinHd,'color','r'); % line색을 red로
set(gcf,'color','w'); % 배경색을 white로
grid on   % 그리드 
% 레벨만들기
xlabel('Salinity (PSU)');
ylabel('Depth (m)');
title('TED6 염분 수직분포도');