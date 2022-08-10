clc;clear all;close all

[f, p]=uigetfile('*.*','select line data');
filedir=[p,f];

data=load(filedir);

depmin=0;depmax=max(data(:,3));
dismin=0;dismax=max(data(:,8));
Xp=[dismin:0.1:dismax];Yp=[depmin:0.1:depmax];
[X,Y]=meshgrid(Xp,Yp);
Zt=griddata(data(:,8),data(:,3),data(:,4),X,Y);
Zs=griddata(data(:,8),data(:,3),data(:,6),X,Y);
Zd=griddata(data(:,8),data(:,3),data(:,7)-1000,X,Y);

leng=length(data);
tick_la_data=[0];
dep_data=[];
for itic=1:1:leng-1
    if data(itic,8)==data(itic+1,8)
        continue
    else
        tick_la_data=[tick_la_data,data(itic+1,8)];
        dep_data=[dep_data,data(itic,3)];
    end
end
dep_data=[dep_data,data(end,3)];
tick_data=tick_la_data./max(tick_la_data);

%--------------------------------------------------------------------------
% Temperature

figure(1);
figure('Position',[10 50 750 350])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7.5 4])
hold on;
[cs,h]=contour(X,-Y,Zt,[20:0.1:22],'linecolor',[.01 .01 .01]);
% contour(X,-Y,Zt,[19.75:0.5:22],'linecolor',[.01 .01 .01]);
pcolor(X,-Y,Zt);colorbar vert;

H1=colorbar('vert');
set(get(H1,'title'),'string','Temp.(℃)');caxis([20 22]);
shading interp;

z_patch=ones(length(dep_data)+2,1);
x_patch=[tick_la_data,max(tick_la_data),0];
y_patch=[-dep_data,-max(dep_data),-max(dep_data)];
patch(x_patch,y_patch,z_patch,[.99 1 1],'edgecolor',[.99 1 1]);

for i_st=1:1:length(dep_data)
    plot(tick_la_data(1,i_st),-[0:5:dep_data(1,i_st)],'marker','square','markeredgecolor','k','markerfacecolor','k','markersize',2);
end

ax(1)=gca;
set(ax(1),'position',[0.07 0.15 0.8 0.73],...
   'Color','none','LineWidth',2,...
   'TickDir','out','box','on');
clabel(cs,h,'fontsize',9,'labelspacing',700,'fontweight','bold');
set(get(ax(1),'xlabel'),'string','Distance(km)');
set(get(ax(1),'ylabel'),'string','Depth(m)');

text(15.5,-30,1,'High tide','fontsize',10,'fontweight','bold');text(15.5,-33,1,'2012.10.12','fontsize',10,'fontweight','bold')
% text(14.5,-29,1,'Low tide','fontsize',10,'fontweight','bold');text(14.5,-32,1,'2012.10.12','fontsize',10,'fontweight','bold')
tit_name=['Gwangwang Temperature'];
title(tit_name,'fontsize',14,'fontweight','bold');
filename=['121012-1 Gwangwang Temperature'];  %% 파일 이름 변경    
print('-dpng',filename);

%--------------------------------------------------------------------------
% Salinity

figure(2);
figure('Position',[10 50 750 350])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7.5 4])
hold on;
% contour(X,-Y,Zs,[29:0.5:31.5],'linecolor',[.01 .01 .01]);
[cs,h2]=contour(X,-Y,Zs,[28.75:0.1:31.5],'linecolor',[.01 .01 .01]);
pcolor(X,-Y,Zs);

H2=colorbar('vert');
set(get(H2,'title'),'string','Sal.(‰)');caxis([29 31.5]);
shading interp;

z_patch=ones(length(dep_data)+2,1);
x_patch=[tick_la_data,max(tick_la_data),0];
y_patch=[-dep_data,-max(dep_data),-max(dep_data)];
patch(x_patch,y_patch,z_patch,[.99 1 1],'edgecolor',[.99 1 1]);

for i_st=1:1:length(dep_data)
    plot(tick_la_data(1,i_st),-[0:5:dep_data(1,i_st)],'marker','square','markeredgecolor','k','markerfacecolor','k','markersize',2);
end

ax(1)=gca;
set(ax(1),'position',[0.07 0.15 0.8 0.73],...
   'Color','none','LineWidth',2,...
   'TickDir','out','box','on');
clabel(cs,h2,'fontsize',9,'labelspacing',700,'fontweight','bold');
set(get(ax(1),'xlabel'),'string','Distance(km)');
set(get(ax(1),'ylabel'),'string','Depth(m)');

text(15.5,-30,1,'High tide','fontsize',10,'fontweight','bold');text(15.5,-33,1,'2012.10.12','fontsize',10,'fontweight','bold')
% text(14.5,-29,1,'Low tide','fontsize',10,'fontweight','bold');text(14.5,-32,1,'2012.10.12','fontsize',10,'fontweight','bold')
tit_name=['Gwangwang Salinity'];
title(tit_name,'fontsize',14,'fontweight','bold');
filename=['121012-1 Gwangwang Salinity'];  %% 파일 이름 변경    
print('-dpng',filename);

%--------------------------------------------------------------------------
% density

figure(3);
figure('Position',[10 50 750 350])
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 7.5 4])
hold on;
% contour(X,-Y,Zd,[20.5:0.1:21.5],'linecolor',[.01 .01 .01]);
[cs,h3]=contour(X,-Y,Zd,[20:0.05:22],'linecolor',[.01 .01 .01]);
pcolor(X,-Y,Zd);

H3=colorbar('vert');
set(get(H3,'title'),'string','σ(kg/m^3)');caxis([20 22]);
shading interp;

z_patch=ones(length(dep_data)+2,1);
x_patch=[tick_la_data,max(tick_la_data),0];
y_patch=[-dep_data,-max(dep_data),-max(dep_data)];
patch(x_patch,y_patch,z_patch,[.99 1 1],'edgecolor',[.99 1 1]);

for i_st=1:1:length(dep_data)
    plot(tick_la_data(1,i_st),-[0:5:dep_data(1,i_st)],'marker','square','markeredgecolor','k','markerfacecolor','k','markersize',2);
end

ax(1)=gca;
set(ax(1),'position',[0.07 0.15 0.8 0.73],...
   'Color','none','LineWidth',2,...
   'TickDir','out','box','on');
clabel(cs,h3,'fontsize',9,'labelspacing',700,'fontweight','bold');
set(get(ax(1),'xlabel'),'string','Distance(km)');
set(get(ax(1),'ylabel'),'string','Depth(m)');

text(15.5,-30,1,'High tide','fontsize',10,'fontweight','bold');text(15.5,-33,1,'2012.10.12','fontsize',10,'fontweight','bold')
% text(14.5,-29,1,'Low tide','fontsize',10,'fontweight','bold');text(14.5,-32,1,'2012.10.12','fontsize',10,'fontweight','bold')
tit_name=['Gwangwang Density'];
title(tit_name,'fontsize',14,'fontweight','bold');
filename=['121012-1 Gwangwang Density'];  %% 파일 이름 변경    
print('-dpng',filename);


%--------------------------------------------------------------------------
% ctd profile
% fig_count=2;
% plot_data=[];
% pdat_num=0;
% xlabels{1} = 'Salinity(‰)';
% xlabels{2} = 'Temperature (℃)';
% ylabels{1} = 'Depth(m)';
% ylabels{2} = 'Depth(m)';
% 
% for ip=1:1:length(data)
%     if data(ip,3)~=0 || ip==1 && ip ~= length(data)
%         pdat_num=pdat_num+1;
%         plot_data(pdat_num,:)=data(ip,:);
%         else if data(ip,3)==0 && ip ~= 1  
%             fig_count=fig_count+1;
%             figure(fig_count);
%             [AX,h1,h2]=plotxx(plot_data(:,6),-plot_data(:,3),plot_data(:,4),-plot_data(:,3),xlabels,ylabels);
%             set(gca,'yGrid','on');
%             set(AX(2),'xlim',([22.8 23.5]));
%             set(AX(1),'xlim',([30.83 31.1]));
%             set(get(AX(1),'xlabel'),'fontweight','bold','fontsize',11)
%             set(get(AX(2),'xlabel'),'fontweight','bold','fontsize',11)
%             profile_name=['CTD profile ',num2str(fig_count-2)];
%             title(profile_name,'fontweight','bold','fontsize',13);
%             print('-dpng',profile_name);
%             plot_data=[];
%             pdat_num=1;
%             plot_data(pdat_num,:)=data(ip,:);
%         end
%     end
%     if ip==length(data)
%     pdat_num=pdat_num+1;
%     plot_data(pdat_num,:)=data(ip,:);        
%     fig_count=fig_count+1;
%     figure(fig_count);
%     [BX,h1,h2]=plotxx(plot_data(:,6),-plot_data(:,3),plot_data(:,4),-plot_data(:,3),xlabels,ylabels);
%     set(BX(2),'xlim',([22.8 23.5]));
%     set(BX(1),'xlim',([30.83 31.1]));
%     set(get(BX(1),'xlabel'),'fontweight','bold','fontsize',11)
%     set(get(BX(2),'xlabel'),'fontweight','bold','fontsize',11)
%     set(gca,'yGrid','on');
%     profile_name=['CTD profile ',num2str(fig_count-2)];
%     title(profile_name,'fontweight','bold','fontsize',13);
%     print('-dpng',profile_name);
% 
%     end
% end

