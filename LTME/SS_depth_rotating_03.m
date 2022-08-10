%% depth 자료로 cutting
for i = 1:length(depth)
    u(i,depth(i)-1:end) = NaN;
    v(i,depth(i)-1:end) = NaN;
    w(i,depth(i)-1:end) = NaN;
    mag(i,depth(i)-1:end) = NaN;      
end

%% shadowzone 판단용
% for i = 1:length(depth)
%     u(i,depth(i)+4:20) = NaN;
%     v(i,depth(i)+4:20) = NaN;
%     mag(i,depth(i)+4:20) = NaN;  
% end

u = flipud(rot90(u(:,1:15)/10));
v = flipud(rot90(v(:,1:15)/10));
w = flipud(rot90(w(:,1:15)/10));
mag = flipud(rot90(mag(:,1:end)/10));

mag(mag == -32768) = nan;
mean_u = nanmean(u);
mean_v = nanmean(v);
mean_w = nanmean(w);

%% rotate angle
angleV = -45;
theta = dir;
%Magnitude remains constant.  Convert to new theta

theta2 = theta - angleV;
%Compensate if new angle has passed 0 degrees
over_ang = find(theta2 > 360);
if isempty(over_ang) == 0;
    theta2(over_ang) = theta2(over_ang)-360;
end

%Reconfigure to new coordinate system
up = mag.*sind(theta2');
across = mag.*cosd(theta2');

% figure()
% mean_up = nanmean(up);
% plot(mean_up);
% hold on;
% plot([0:30:4531],0,'k.','markersize',2)
% set(gca,'xTick',[0:144:4752]);                
% set(gca,'xTickLabel',[0:1:33]);
% xlim([0 4531]); %2014년 관측과 기간을 맞춰봄
% ylim([-150 100]);

figure('name','check:ratate angle');
subplot(211);
contourf(across); shading flat;
ylabel('Across','fontsize',14);
% hold on; plot(dd,'linewidth',2)
caxis([-50 50]);colorbar
set(gca,'xTick',[0:144:8928],'xTickLabel',[6:31,1:10],'xlim',[0 17*144]);    
set(gca,'yTick',[-1:2:15],'yTickLabel',[0:2:20],'ylim',[-1 14]); 

subplot(212);
contourf(up); shading flat;
ylabel('Up','fontsize',14);
caxis([-80 80]);colorbar
set(gca,'xTick',[0:144:8928],'xTickLabel',[6:31,1:10],'xlim',[0 17*144]);    
set(gca,'yTick',[-1:2:15],'yTickLabel',[0:2:20],'ylim',[-1 14]); 
% hold on; plot(depth,'linewidth',2)