sp_mean_no3_yr
sp_mean_temp_yr
sp_mean_salt_yr
%% salt
clearvars yp_w_salt
color_pic = lines(size(regime_salt_yr,1));
marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<',...
    'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<'};
xp = 1:22;
j=0
figure; hold on;
for i = 1:size(sp_mean_salt_yr,1)
  clearvars reg_data_salt xp_w_salt pf_w_salt
reg_data_salt = sp_mean_salt_yr(i,:);
if isnan(reg_data_salt(1)) == 0 
    j = j+1;
    xp_w_salt = find(isnan(reg_data_salt)==0);
    pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
    yp_w_salt(i,:) = polyval(pf_w_salt,xp);
    scatter(1:22,sp_mean_salt_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_salt(i,:),'color',color_pic(i,:));
    sp_coeff_salt(i,:) = pf_w_salt;
    sp_temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM-Ç¥Ãþ¿°ºÐ ¿¬Æò±Õ(°ø°£Æò±Õ)','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])


%% no3
clearvars yp_w_no3
j=0
figure; hold on;
for i = 1:size(sp_mean_no3_yr,1)
  clearvars reg_data_no3 xp_w_no3 pf_w_no3
reg_data_no3 = sp_mean_no3_yr(i,:);
if isnan(reg_data_no3(1)) == 0 
    j = j+1;
    xp_w_no3 = find(isnan(reg_data_no3)==0);
    pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
    yp_w_no3(i,:) = polyval(pf_w_no3,xp);
    scatter(1:22,sp_mean_no3_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_no3(i,:),'color',color_pic(i,:));
    sp_coeff_no3(i,:) = pf_w_no3;
    sp_temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM°üÃø-Ç¥ÃþÁú»ê¿° ¿¬Æò±Õ(°ø°£Æò±Õ)','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%temp
clearvars yp_w_temp
j=0
figure; hold on;
for i = 1:size(sp_mean_temp_yr,1)
  clearvars reg_data_temp xp_w_temp pf_w_temp
reg_data_temp = sp_mean_temp_yr(i,:);
if isnan(reg_data_temp(1)) == 0 
    j = j+1;
    xp_w_temp = find(isnan(reg_data_temp)==0);
    pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
    yp_w_temp(i,:) = polyval(pf_w_temp,xp);
    scatter(1:22,sp_mean_temp_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_temp(i,:),'color',color_pic(i,:));
    sp_coeff_temp(i,:) = pf_w_temp;
    sp_temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('temperature (^oC)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM°üÃø-Ç¥Ãþ¼ö¿Â ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
ylim([13 20])

