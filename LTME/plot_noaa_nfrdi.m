clc; clear all; close all;

data = load('nfrdi_trend_with_noaa.txt');

x = data(:,4); % 수과원 자료
y = data(:,5); % noaa 자료

figure('Position', [100, 100, 400, 400])
plot(x,y,'ko','markersize',5,'markerfacecolor','k');
xlabel('NFRDI','fontsize',14);
ylabel('NOAA/AVHRR','fontsize',14);
axis([-0.02 0.04 -0.02 0.04]);
hold on;
plot(-0.02:0.01:0.04, -0.02:0.01:0.04,'k-')
axis square

%--- save figure --------------------------------------
%     outpath = '.\figure\';
%     out=[outpath,'NFRDI_vs_KOHA'];
%     set(gcf,'renderer','painter');
%     set(gcf, 'PaperUnits', 'inches');
%     x_width=4;
%     y_width=4;
%     set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
%     saveas(gcf,out,'tif');
%------------------------------------------------------
