

%% salt
color_pic = {'b','g','r','k','m'}
xp = 1:40;
figure; hold on;
for i = 1:5
  clearvars reg_data_salt xp_w_salt pf_w_salt
reg_data_salt = reg_clim_salt(i,:);
xp_w_salt = find(isnan(reg_data_salt)==0);
pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
yp_w_salt(i,:) = polyval(pf_w_salt,xp);
scatter(1:40,reg_clim_salt(i,:),color_pic{i});
plot(1:40, yp_w_salt(i,:),color_pic{i});
coeff_salt(i,:) = pf_w_salt;
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥Ãþ¿°ºÐ ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([32 35])


%% no3
color_pic = {'b','g','r','k','m'}
figure; hold on;
for i = 1:5
clearvars reg_data_no3 xp_w_no3 pf_w_no3
if i ~= 2 && i ~= 4 
reg_data_no3 = reg_clim_no3(i,:);
xp_w_no3 = find(isnan(reg_data_no3)==0);
pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
yp_w_no3(i,:) = polyval(pf_w_no3,xp);
scatter(1:40,reg_clim_no3(i,:),color_pic{i});
plot(1:40, yp_w_no3(i,:),color_pic{i});
coeff_no3(i,:) = pf_w_no3;
hold on
end
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥ÃþÁú»ê¿° ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%temp
color_pic = {'b','g','r','k','m'}
xp = 1:40;
figure; hold on;
for i = 1:5
  clearvars reg_data_temp xp_w_temp pf_w_temp
reg_data_temp = reg_clim_temp(i,:);
xp_w_temp = find(isnan(reg_data_temp)==0);
pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
yp_w_temp(i,:) = polyval(pf_w_temp,xp);
scatter(1:40,reg_clim_temp(i,:),color_pic{i});
plot(1:40, yp_w_temp(i,:),color_pic{i});
coeff_temp(i,:) = pf_w_temp;
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:5:40]);
set(gca,'xlim',[1 40]);
set(gca,'xticklabel',1980:5:2019);
title('Á¤¼±°üÃø-Ç¥Ãþ¼ö¿Â ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)


