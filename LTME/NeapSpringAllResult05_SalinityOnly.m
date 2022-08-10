%% salinity vertical profile
st_st = 6;
figure('position', [10, 10, 600, 1100])
subplot(411)
data = neap_ebb_sal(2:end,:);
[a b]=size(data);
for i = 1:b-1
    l = sum(~isnan(data(:,i)));
%     plot(data(1:l,i+1),data(1:l),'linewidth',2,'color',[255-8*i 10+7*i 12*i]/255)
    plot(data(:,i),'linewidth',2,'color',[255-8*(b-i) 8*(b-i) 255-8*i]/255)
    hold on;
end
 set(gca,'ylim',[0 35],'xlim',[0 40]);
 xlabel('Depth (m)','fontsize',14);
 ylabel('Salinity (¢¶)','fontsize',14);
%  legend('st1','st2','st3','st4','st5','st6','st7','st8','st9','st10','st11','st12','location','eastoutside')
view([90 90]);
text(33, 7, 'Neap Ebb','fontsize',14);

subplot(412)
data = neap_flood_sal(2:end,:);
[a b]=size(data);
for i = 1:b-1
    l = sum(~isnan(data(:,i)));
%     plot(data(1:l,i+1),data(1:l),'linewidth',2,'color',[255-8*i 10+7*i 12*i]/255)
    plot(data(:,i),'linewidth',2,'color',[255-8*(b-i) 8*(b-i) 255-8*i]/255)
    hold on;
end
 set(gca,'ylim',[0 35],'xlim',[0 40]);
 xlabel('Depth (m)','fontsize',14);
 ylabel('Salinity (¢¶)','fontsize',14);
%  legend('st1','st2','st3','st4','st5','st6','st7','st8','st9','st10','st11','st12','location','eastoutside')
view([90 90]);
text(33, 7, 'Neap Flood','fontsize',14);

subplot(413)
data = spring_ebb_sal(2:end,:);
[a b]=size(data);
for i = 1:b-1
    l = sum(~isnan(data(:,i)));
%     plot(data(1:l,i+1),data(1:l),'linewidth',2,'color',[255-8*i 10+7*i 12*i]/255)
    plot(data(:,i),'linewidth',2,'color',[255-8*(b-i) 8*(b-i) 255-8*i]/255)
    hold on;
end
 set(gca,'ylim',[0 35],'xlim',[0 40]);
 xlabel('Depth (m)','fontsize',14);
 ylabel('Salinity (¢¶)','fontsize',14);
%  legend('st1','st2','st3','st4','st5','st6','st7','st8','st9','st10','st11','st12','location','eastoutside')
view([90 90]);
text(33, 7, 'Spring Ebb','fontsize',14);

subplot(414)
data = spring_flood_sal(2:end,:);
[a b]=size(data);
for i = 1:b-1
    l = sum(~isnan(data(:,i)));
%     plot(data(1:l,i+1),data(1:l),'linewidth',2,'color',[255-8*i 10+7*i 12*i]/255)
    plot(data(:,i),'linewidth',2,'color',[255-8*(b-i) 8*(b-i) 255-8*i]/255)
    hold on;
end
 set(gca,'ylim',[0 35],'xlim',[0 40]);
 xlabel('Depth (m)','fontsize',14);
 ylabel('Salinity (¢¶)','fontsize',14);
%  legend('st1','st2','st3','st4','st5','st6','st7','st8','st9','st10','st11','st12','location','eastoutside')
view([90 90]);
text(33, 7, 'Spring Flood','fontsize',14);



