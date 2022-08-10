% tbase와 고해상도 중남미 지도  및 정점표현
% tbase_exam3.m
%
clc;clear all;close all;
m_proj('lambert','lon',[-100 -65],'lat',[-20 50]);
m_tbase('contourf',[-5000:500:0],'edgecolor','b');
m_grid('linestyle','none','tickdir','out','linewidth',3);
m_gshhs_h('patch',[.7 .7 .7]);
% m_coast('patch',[.7 .7 .7],'edgecolor','r'); 
% low resolution coast

m_line(-84.17,9.40,'marker','square','color','r');
m_line(-79.56,8.97,'marker','square','color','r');
m_line(-79.90,9.35,'marker','square','color','r');
