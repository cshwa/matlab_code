clc;clear all;close all;

f_data=[];

    [f, p]=uigetfile('*.*','select ADCP data file');
        filedir=[p,f];
        disp(filedir);
    data=xlsread(filedir);data=data(5:end,11:end);
    uvdata=data(:,2:end)/10;
    leng=length(data(:,1));
    
    x=[1:1:15];y=[1:1:leng];
    u=[];v=[];
    for i=1:1:15
    u=[u,uvdata(:,2*i-1)];
    v=[v,uvdata(:,2*i)];
    end
    
%     figure(1);
%     scrsz = get(0,'ScreenSize');
%     figure('Position',[1 scrsz(4)/3.7 scrsz(3)/2 scrsz(4)/1.5])
%  set(gca,'position',[0.05 0.12 0.8 0.5])
%     quiver(x,y,u,v)
%     set(gca,'view',[90 -90])

    
    figure(2)
    hold on;
     x2=[1:0.1:15];y2=[1:0.1:leng];
     [X,Y]=meshgrid(x2,y2);
     
     Zu=griddata(x,y,u,X,Y);
    [cs,h]=contour(X,Y,Zu,[-80,-40,40],'linecolor',[0 0 0]);
    [cs1,h1]=contour(X,Y,Zu,[0 0],'linecolor',[.99 .99 1],'linewidth',1.5);
    pcolor(X,Y,Zu);shading interp;
    H1=colorbar('vert');
    set(get(H1,'title'),'string','Vel.(cm/s)');caxis([-80 40]);

    ax(1)=gca;
set(ax(1),'position',[0.05 0.12 0.8 0.8],'view',[90 -90],...
   'Color','none','LineWidth',2,...
   'TickDir','out','box','on',...
    'ytick',[6:30:156],'yticklabel',{'10/11 14h','19h','10/12 00h','5h','10h','15h'},...
    'xtick',[1:2:15],'xticklabel',{'7','6','5','4','3','2','1','0'});
clabel(cs,h,'fontsize',9,'labelspacing',700,'fontweight','bold');
clabel(cs1,h1,'fontsize',9,'labelspacing',700,'fontweight','bold');
set(get(ax(1),'ylabel'),'string','Time');
set(get(ax(1),'xlabel'),'string','Depth(m)');

tit_name=['ADCP U velocity'];
title(tit_name,'fontsize',14,'fontweight','bold');
filename=['ADCP U velocity'];  %% 파일 이름 변경    
print('-dpng',filename);

    figure(3)
     hold on;
     Zv=griddata(x,y,v,X,Y);
    [cs,h]=contour(X,Y,Zv,[-80,-40,40],'linecolor',[0 0 0]);
    [cs1,h1]=contour(X,Y,Zv,[0 0],'linecolor',[.99 .99 .99],'linewidth',1.5);
    pcolor(X,Y,Zv);shading interp;
    H1=colorbar('vert');
    set(get(H1,'title'),'string','Vel.(cm/s)');caxis([-80 40]);

    ax(1)=gca;
set(ax(1),'position',[0.05 0.12 0.8 0.8],'view',[90 -90],...
   'Color','none','LineWidth',2,...
   'TickDir','out','box','on',...
    'ytick',[6:30:156],'yticklabel',{'10/11 14h','19h','10/12 00h','5h','10h','15h'},...
     'xtick',[1:2:15],'xticklabel',{'7','6','5','4','3','2','1','0'});
clabel(cs,h,'fontsize',9,'labelspacing',700,'fontweight','bold');
clabel(cs1,h1,'fontsize',9,'labelspacing',700,'fontweight','bold');

set(get(ax(1),'ylabel'),'string','Time');
set(get(ax(1),'xlabel'),'string','Depth(m)');

tit_name=['ADCP V velocity'];
title(tit_name,'fontsize',14,'fontweight','bold');
filename=['ADCP V velocity'];  %% 파일 이름 변경    
print('-dpng',filename);