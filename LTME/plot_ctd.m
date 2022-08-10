clc;clear all;close all;

data=load('timeseris_ctd.dat');

max_dep=max(data(:,2));

time_num=26;
n=0;
for i=1:1:length(data)-1
    if data(i,1)~=data(i+1,1) || i==length(data)-1
        max2_dep=max(data(i-n:i,2));
        gap_dep=max_dep-max2_dep;
        data(i-n:i,2)=data(i-n:i,2)+gap_dep;
        n=0;
    else
        n=n+1;
    end
end

dep_cutdata=[data(1,2)];
for itic=1:1:length(data)-1
    if data(itic,1)==data(itic+1,1)
        continue
    else
        dep_cutdata=[dep_cutdata,data(itic+1,2)];
    end
end
     
xmin=min(data(:,1));xmax=max(data(:,1));
ymin=min(data(:,2));ymax=max(data(:,2));
X=[xmin:0.05:xmax];Y=[ymin:0.05:ymax];
[Xi,Yi]=meshgrid(X,Y);
Zt=griddata(data(:,1),data(:,2),data(:,3),Xi,Yi);
Zs=griddata(data(:,1),data(:,2),data(:,4),Xi,Yi);

%--------------------------------------------------------------------------
%temp

figure(1);
hold on;
[cs,h]=contour(Xi,-Yi,Zt,[18:0.5:23],'linecolor',[0 0 0]);
pcolor(Xi,-Yi,Zt);colorbar vert;

H1=colorbar('vert');
set(get(H1,'title'),'string','Temp.(℃)');caxis([18 23]);
shading interp;

z_patch=ones(length(dep_cutdata)+2,1);
x_patch=[1.,1:1:time_num,time_num];
y_patch=[0.1,dep_cutdata,0.1]*(-1);
patch(x_patch,y_patch,z_patch,[0.99 1 1],'edgecolor',[0.99 1 1]);
ylim([max(data(:,2))*(-1) 0]);


ax(1)=gca;
set(ax(1),'position',[0.12 0.12 0.75 0.70],...
   'Color','none','LineWidth',2,...
   'TickDir','out','box','on',...
    'xtick',[1:5:26],'xticklabel',{'10/11 14h','19h','10/12 00h','5h','10h','15h'},...
    'xminortick','on');
clabel(cs,h,'fontsize',9,'labelspacing',700,'fontweight','bold');

set(get(ax(1),'xlabel'),'string','Time');
set(get(ax(1),'ylabel'),'string','Depth(m)');

tit_name=['Temperature variation'];
title(tit_name,'fontsize',14,'fontweight','bold');
filename=['GwangYang variation temp '];  %% 파일 이름 변경    
print('-dpng',filename);

%--------------------------------------------------------------------------
%salinity

figure(2);
hold on;
[cs1,h1]=contour(Xi,-Yi,Zs,[12:1:28],'linecolor',[0 0 0]);
pcolor(Xi,-Yi,Zs);colorbar vert;

H1=colorbar('vert');
set(get(H1,'title'),'string','Sal.(psu)');caxis([12 28]);
shading interp;

z_patch=ones(length(dep_cutdata)+2,1);
x_patch=[1.,1:1:time_num,time_num];
y_patch=[0.1,dep_cutdata,0.1]*(-1);
patch(x_patch,y_patch,z_patch,[0.99 1 1],'edgecolor',[0.99 1 1]);
ylim([max(data(:,2))*(-1) 0]);



ax(1)=gca;
set(ax(1),'position',[0.12 0.12 0.75 0.70],...
   'Color','none','LineWidth',2,...
   'TickDir','out','box','on','Xminortick','on',...
   'xtick',[1:5:26],'xticklabel',{'10/11 14h','19h','10/12 00h','5h','10h','15h'},...
   'xminortick','on');
clabel(cs1,h1,'fontsize',9,'labelspacing',800,'fontweight','bold');

set(get(ax(1),'xlabel'),'string','Time');
set(get(ax(1),'ylabel'),'string','Depth(m)');

tit_name=['Salinity variation'];
title(tit_name,'fontsize',14,'fontweight','bold');
filename=['GwangYang variation sal '];  %% 파일 이름 변경    
print('-dpng',filename);

