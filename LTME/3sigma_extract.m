clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f w_chl_mi2
% raw

idx1 = find(isnan(w_chl_mi) == 0);
w_chl_mi2=w_chl_mi;
w_chl_mi2(find(w_chl_mi > 3*std(w_chl_mi(idx1))))=NaN;
 figure; hold on;
 plot(w_chl,'.');
 xlabel('time(year)','fontsize',13)
ylabel('chl (mg/L)','fontsize',13)
title(['sumjin(songjung)-polynomial chl raw - climate'],'fontsize',13)
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[1 t_tick(4:4:end)]);
set(gca,'xlim',[1 t_tick(end)]);
set(gca,'xticklabel',1980:4:2020);
%regression
%slope y = b1*x
nonan=w_chl_mi2(~isnan(w_chl_mi2));
idx = find(isnan(w_chl_mi2) == 0);
nanaxisx=find(isnan(w_chl_mi2) == 1);
b1 = idx'\nonan'; % x\y for getting slop
yCalc1 = b1*idx;
% Slope & Intercept y = b0 + b1*x
X = [ones(length(idx'),1) idx']; %b_0 b_1 
b = X\nonan';
yCalc2 = [1:length(w_chl_mi)].*b(2) + b(1);
xp = [1:length(w_chl_mi)];
pf_w_chl = polyfit(idx,w_chl_mi2(idx),3);
poly_w_chl = polyval(pf_w_chl,xp);
for i = 1:length(yCalc2)
        w_chl_recon(i) = yCalc2(i)  + yp_w_chl(indx_366{i});
        w_chl_recon_poly(i) = poly_w_chl(i)  + yp_w_chl(indx_366{i});
end
plot([1:length(w_chl_mi)], w_chl_recon,'--','color','r','linew',1)
plot([1:length(w_chl_mi)], poly_w_chl,'--','color','c','linew',1)
plot([1:length(w_chl_mi)], w_chl_recon_poly,'--','color','g','linew',1)
ylim([0 5*std(w_chl_mi(idx1))])
plot([1:length(w_chl_mi)],yCalc2,'--','color','m','linew',2)
gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
% print('-dpng', ['sumjin(songjung)_chl_raw&recon']);
set(gca,'xlim',[1 t_tick(10)]);
set(gca,'xtick',[1 t_tick(1:1:end)]);
set(gca,'xticklabel',1980:1:2020);
% print('-dpng', ['sumjin(songjung)_chl_raw&recon-80s']);