% tbase를 이용한 우리나라 부근 수심도
% tbase_exam2.m
clc;clear all;close all;
m_proj('mercator','lon',[120 140],'lat',[30 50]);
m_tbase('contourf');
% 위경도 표시방법
m_grid('linestyle','none','tickdir','in','linewidth',2);
m_grid('box','fancy','tickdir','in');