clc;close all;clear all;

[f, p]=uigetfile('*.*','select station coordinate');
filedir=[p,f];

fid=fopen(filedir);
c_header=textscan(fid,'%s %s %s %s',1);
c_data=textscan(fid,'%s %s %f %f');
fclose(fid);
data=[c_data{3},c_data{4}];

m_proj('mercator','long',[127.7 128],'lat',[34.9 35.05]);
m_gshhs_h('color','k');
m_gshhs_h('patch',[.66 .66 .66]);
m_grid('box','fancy','tickdir','in','linewidth',1);
m_text(127.87,34.92,'Namhae','fontsize',15);

% m_line(127.7790,34.9908,'marker','square','markeredgecolor','r','markerfacecolor','r','markersize',3);
% m_text(127.7790,34.9908+0.01,'ADCP','fontsize',10,'fontweight','bold','color','r');
% 
% 
% m_line(127.7774,34.9884,'marker','square','markeredgecolor','k','markerfacecolor','k','markersize',3);
% m_text(127.7774,34.9884+0.01,'CTD','fontsize',10,'fontweight','bold','color','k');

m_line(data(1:8,2),data(1:8,1),'color','r','linewi',1.5,'marker','square','markeredgecolor','r','markerfacecolor','r','markersize',3);
for i=1:1:8
st_n=['high',num2str(i)];
m_text(data(i,2),data(i,1)+0.01,st_n,'fontsize',10,'fontweight','bold','color','r');
end

m_line(data(9:16,2),data(9:16,1),'color','k','linewi',1.5,'marker','square','markeredgecolor','k','markerfacecolor','k','markersize',3);

for i=9:1:16
st_n=['low',num2str(i-8)];
m_text(data(i,2),data(i,1)-0.01,st_n,'fontsize',10,'fontweight','bold');
end

print('-dpng','ctd station')