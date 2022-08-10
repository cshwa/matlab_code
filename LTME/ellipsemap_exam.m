% 조류타원도 그리기
% m-file: ellipsemap_exam.m
% using m_ellipse.m
% 다음은 조류조화분해 결과
% NAME      SPEED       MAJOR   MINOR   INC    G      G+     G-
% 13 M2    0.08051140     0.730   0.088   68.9   61.5  352.6  130.4  00ic-1.dat
%  6 M2    0.08051140     0.655  -0.009   41.2   50.4    9.1   91.6  02ic-2.dat

close all;clear all;clc
% 지도작성
m_proj('mercator','lat',[37.0 37.6],'long',[126.15 126.9]);
m_gshhs_f('color','k');hold on    % 초고해상도 지도
m_gshhs_f('patch',[.8 .8 .8]);    % 지도색 조정
m_grid('linestyle','none','box','fancy','tickdir','in');  % 위도 및 경도선 그리기
hold on
% m_ellipse를 활용한 tidal ellipse 작성
scale = 15;                                    %  기준 Scale조정
maj1=0.730/scale;min1=0.088/scale;inc1=68.9;    %  장축
inc1=inc1*pi/180;
maj2=0.655/scale;min2=-0.009/scale;inc2=41.2;   % 단축
inc2=inc2*pi/180;
m_ellipse(maj1,min1,inc1,126.4097,37.2781,'r'); % 장축,단축,각도,중심좌표
hold on
m_ellipse(maj2,min2,inc2,126.5222,37.3736,'r');  % 장축,단축,각도,중심좌표
set(gcf,'Color','w')