clear all; clc; close all;


% list = dir(f);

transp = zeros(10,6);

for i = 1:10
f = strcat('roms_transport_monthly_csh_',num2str(i),'.txt');

temp= load(f);
f
transp(i,:) = temp;

end

save('RCP_2yr_2007to2100_v2.mat','transp');




figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;
out_name = ['all_section_transport_RCP85'];

figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)  
hold on;
plot( transp(3:end,1),'linewidth',2,'color','g');
plot( transp(3:end,2),'linewidth',2,'color','b');
plot( transp(3:end,3),'linewidth',2,'color','k');
plot( transp(3:end,4),'linewidth',2,'color','c');
plot( transp(3:end,5),'linewidth',2,'color','r');
plot( transp(3:end,6),'linewidth',2,'color','m');
% legend('r-mouth(low)', 'r-mouth(high)','l-mouth', 'left(ent)', 'main-channel', 'jinju-channel');
%         ylim([34.8 34.96])
xlim([1 8]);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
xlabel('year','color','k','FontSize',fontSizeLab,'fontweight','bold')
ylabel('Transport(CMS)','color','k','FontSize',fontSizeLab,'fontweight','bold')
set(gca, 'XTick',(1:1:10))
set(gca, 'XTicklabel',2021:10:2091)
print(gcf,['figures\',out_name,'_vertical_integral','.png'],'-dpng','-r200');


%%%

nname={'r-mouth(low)', 'r-mouth(high)','l-mouth', 'left(ent)', 'main-channel', 'jinju-channel'};
ccolor_p = {'g','b','k','c','r','m'}; 
figPos = [0 0 5 4];
fontSizeTick = 12;
fontSizeLab  = 12;
fontSizeCLab = 12;

for i = 1:length(nname)
out_name = [nname{i},'_section_transport_RCP85'];

figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos)  
hold on;
plot( transp(3:end,i),'linewidth',2,'color',ccolor_p{i});
% legend('r-mouth(low)', 'r-mouth(high)','l-mouth', 'left(ent)', 'main-channel', 'jinju-channel');
%         ylim([34.8 34.96])
xlim([1 8]);
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on');
xlabel('year','color','k','FontSize',fontSizeLab,'fontweight','bold')
ylabel('Transport(CMS)','color','k','FontSize',fontSizeLab,'fontweight','bold')
set(gca, 'XTick',(1:1:10))
set(gca, 'XTicklabel',2021:10:2091)
print(gcf,['figures\',out_name,'_vertical_integral','.png'],'-dpng','-r200');
hold off;
end
