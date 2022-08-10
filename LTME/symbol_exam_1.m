% 지도에 심볼표시
% symbol_exam.m
%
clc;clear all;close all;
m_proj('mercator','lon',[128 134],'lat',[33 39]);  % 경도128~134, 위도 33~39지역의도를 작성
m_gshhs_h('patch',[1 1 .7]);                       % 지도 배경색상 조정
m_grid('linestyle','none','tickdir','in','linewidth',1)  % 위경도를 표시
m_text(130.5,38.5,'EAST SEA','vertical','top');          % 지명을 지도위에 표시
m_text(131.8,34.55,'JAPAN','vertical','top');m_text(128.8,37.5,'Donghae','vertical','top');
m_text(129.2,36.5,'Hupo','vertical','top');m_text(129.06,36.1,'Pohang','vertical','top');
m_text(129.2,35.6,'Ulsan','vertical','top');m_text(129,35.25,'Busan','vertical','top');
m_text(131.85,37.43,'Dokdo','vertical','top')
m_line(131,36,'marker','square','markersize',8,'color','r');     % 빨간네모 심볼
m_line(130.8,37,'marker','diamond','markersize',8,'color','g');  % 녹색 마름모 심볼