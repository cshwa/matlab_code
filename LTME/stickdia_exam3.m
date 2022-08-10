% m-file: stickdia_exam3.m
% By peter
% 주의: x축에 Date number표현이 미비함

clc;clear all;close all;
a = textread('file12_Nc.dat','');
Data1 = a(:,6:17);              % 원하는 부분만 추출 후 행렬작성
Nc = reshape(Data1',[],1)       % 행렬DATA를 1열 행렬DATA로 변환

% East Component Extract and make matrix
b = textread('file13_Ec.dat','');
Data2 = b(:,6:17);              % 원하는 부분만 추출 후 행렬작성
Ec = reshape(Data2',[],1)       % Data2' 는 Transpose를 의미, 1열 자료로 변환

% date number 
n=length(Ec);  % 자료의 갯수
dn=1:n;
dn=dn';  % 1열 Data로 변경
FData = [dn,Nc,Ec]                 %Fdata 행렬표시
fid = fopen('Finaldata.dat','w');
fprintf(fid,'%6d%6d%6d\n',FData');   % FData 형식 그대로 출력
fclose(fid)
hold on                 % 꼭 기록해야 함
% plot 함수를 이용한 Stick diagram 작성
for i = 1:n;
x1 = dn(i);y1 = 0;
x2 = dn(i)+Ec(i);y2 = Nc(i);
plot([x1 x2],[y1 y2],'b');hold on;
end
hold off
% Label 이름지정하기 
title('Stick Diagram using plot function');
xlabel(' Time number');
ylabel('Speed(cm/sec)')
axis equal
axis([0 700 -100 100]);    % AXIS([XMIN XMAX YMIN YMAX])