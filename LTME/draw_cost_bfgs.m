load map_costf.mat;
load bfgs_ana.mat;

for i=1:29
    xn(i) = bfgs_ana(i).x(1);
    zn(i) = bfgs_ana(i).x(3);
    pn(i) = bfgs_ana(i).x(4);
    costf(i) = bfgs_ana(i).cost;
end

figure(1);
contour(xi,zi,cost,[0:100:20000],'fill','on');
colorbar;
caxis([0 20000]);
xlabel('x','fontsize',15);
ylabel('z','fontsize',15);
set(gcf,'position',[0 0 500 400]);
hold on;
plot(1,5,'r+');
plot(xn,zn,'m+');

figure(2);
plot(costf);
set(gcf,'position',[500 0 500 400]);
xlabel('assimilation step','fontsize',15);
ylabel('cost function','fontsize',15);

figure(3);
plot(pn);
set(gcf,'position',[500 200 500 400]);
xlabel('assimilation step','fontsize',15);
ylabel('Parameter P','fontsize',15);
