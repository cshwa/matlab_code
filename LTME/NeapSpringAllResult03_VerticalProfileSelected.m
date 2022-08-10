%% 입구에서 상류 비교
en = 15;
figure()
subplot(341)
contourf(flipud(neap_ebb_temp(2:en,1:10))); shading flat;
caxis([6.5 9.5]);
title('Neap ebb','fontsize',14);
subplot(345)
contourf(flipud(neap_ebb_sal(2:en,1:10))); shading flat;
caxis([0 30]);
subplot(349)
contourf(flipud(neap_ebb_den(2:en,1:10))); shading flat;
caxis([1000 1027]);

subplot(342)
contourf(flipud(neap_flood_temp(2:en,1:10))); shading flat;
caxis([6.5 9.5]);
title('Neap flood','fontsize',14);
subplot(346)
contourf(flipud(neap_flood_sal(2:en,1:10))); shading flat;
caxis([0 30]);
subplot(3,4,10)
contourf(flipud(neap_flood_den(2:en,1:10))); shading flat;
caxis([1000 1027]);


subplot(343)
contourf(flipud(spring_ebb_temp(2:en,1:10))); shading flat;
caxis([6.5 9.5]);
title('Spring ebb','fontsize',14);
subplot(347)
contourf(flipud(spring_ebb_sal(2:en,1:10))); shading flat;
caxis([0 30]);
subplot(3,4,11)
contourf(flipud(spring_ebb_den(2:en,1:10))); shading flat;
caxis([1000 1027]);

subplot(344)
contourf(flipud(spring_flood_temp(2:en,1:10))); shading flat;
caxis([6.5 9.5]); 
title('Neap flood','fontsize',14);
subplot(348)
contourf(flipud(spring_flood_sal(2:en,1:10))); shading flat;
caxis([0 30]);
subplot(3,4,12)
contourf(flipud(spring_flood_den(2:en,1:10))); shading flat;
caxis([1000 1027]);