% Auto Spectrum
% m-file: graphic_exam.m
% Kim Dae Hyun ? 2002. 4. 19.


[filename, procDir] = uigetfile('*.*','Select the Auto spectrum data');
fname = strcat(procDir,filename);
data1= textread(fname,'');
fprintf('Processing %s\n',fname);
   
x_xx=data1(:,1);
y_yy=data1(:,4);
   
index_x = find(x_xx <0.5);
x_x = x_xx(index_x);
y_y = y_yy(index_x);

x=([0.2,0.2]);
y=([620 1960]);

loglog(x,y,'k');hold on;
text(x(1)+0.03,sum(y)/2,'90 %');
[ax,hlT,hlS] = plotyy(x_x,y_y,x_x,y_y);
set(ax(1),'Xscale','log','Yscale','log');
set(ax(1),'XAxisLocation','bottom','YAxisLocation','left');
set(get(ax(1),'XLabel'),'String','Frequency (cph)');
set(get(ax(1),'YLabel'),'String','Spectral Density ((cm/sec)^2/cph)');
set(ax(1),'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
set(ax(1),'YMinorTick','on');
set(hlT,'color',[0 0 0]);

xtick_l=([240,120,48,24,12,6,2]);
x_hour=([1/240,1/120,1/48,1/24,1/12,1/6,1/2]);
set(ax(2),'xtick',x_hour,'Xticklabel',xtick_l);
set(ax(2),'Xscale','log','Yscale','log');
set(ax(2),'XAxisLocation','top','YAxisLocation','right');
set(get(ax(2),'XLabel'),'String','Period (hours)');
set(ax(2),'XMinorTick','off','YTick',[]);
set(ax(2),'Xcolor',[0 0 0],'Ycolor',[0 0 0]);
set(hlS,'color',[0 0 0]);

set(ax,'TickDir','out');
set(ax,'LineWidth',[1.0]);
set(ax,'TickLength',[0.02 0.04]);
set(ax,'Position',[0.15 0.15 0.7 0.7]); 

axis(ax(1),[0.0004 0.5 0.1 1000000]);
axis(ax(2),[0.0004 0.5 0.1 1000000]);
axis(ax(1),'square');
axis(ax(2),'square');

clear;
hold off;
set(gcf,'PaperPosition',[5 9 11 11]);