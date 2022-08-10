cln
set(gcf, 'units', 'norm');
set(gcf, 'position', [.1 .1 .8 .8]);

m_proj('mercator','lon',[128 132],'lat',[35 39]);
load es_skku.mat;
load land1.mat;
mymap = [jet(64);land];

data=xlsread('eastdrilloct');
y=data(:,2);x=data(:,3);t=data(:,4),s=data(:,5);
yg=[36.0:0.005:38.3];
xg=[128.3:0.005:131.5];
[xl,yl]=meshgrid(xg,yg);
temp=griddata(x,y,t,xl,yl,'cubic');

t1 = (temp - 25)./15; 
tp = ([15:0.5:25] - 25)./15;
m_contour(xl,yl,t1,tp,'k');

hold on

z1 = dep./3000;
[a b] = size(z1);
z1c = reshape(z1, a*b, 1);
z2 = reshape(z1c, a, b);
m_pcolor(lon,lat,z2); shading flat;

colormap(mymap);
caxis([-1 1]);
clear data x y t s yg xg temp dep z1 ic z1c lon lat

m_grid('xtick',8,'ytick',9,'fontsize',8,'fontname','times','box','fancy','linewidth',2,'tickdir','in')
xlabel('Longitude','fontname','times','fontsize',12,'fontweight','bold')
ylabel('Latitude','fontname','times','fontsize',12,'fontweight','bold')
print(gcf,'-dtiff','es_basemap.tif')
