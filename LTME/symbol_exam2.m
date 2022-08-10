% 지도상에 grid라인,symbol 및 위경도 위치 변경
% symbol_exam2.m
%
clc;clear all;close all;
m_proj('UTM','lon',[124 130],'lat',[31.5 35]);
m_gshhs_i('color','k'); %해안선 검정
m_gshhs_i('patch',[.7 .7 .7]); %육지부분 흐린흑색지정
% 격자표시 축 조정
m_grid('box','on','tickdir','in','xaxisloc','bottom','yaxisloc','right','clip','on');
m_line(128,33,'marker','square','MarkerEdgeColor','r',...
'MarkerFaceColor','g','MarkerSize',9,'color','r');
m_line(128,33.5,'marker','p','MarkerEdgeColor','k',...
'MarkerFaceColor','r','MarkerSize',10,'color','r');
m_line(128.7,34.15,'marker','o','MarkerEdgeColor','y',...
'MarkerFaceColor','m','MarkerSize',12,'color','r');
set(gcf,'Color','w');              % 배경색 White 