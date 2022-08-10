%% load and plot transport
clearvars utran_mod vtran_mod tran_v_* tran_u_*
for iy = 1982:2020
    clearvars utran_mod vtran_mod 
sav_path = ['J:\fennel_NWP\Output\noep_f\',num2str(iy),'\'];
utran_mod=ncread([sav_path,'roms_monthly_avg.nc'],'Huon');
vtran_mod=ncread([sav_path,'roms_monthly_avg.nc'],'Hvom');

if size(vtran_mod,4) > 12
    utran_mod(:,:,:,13) = [];
    vtran_mod(:,:,:,13) = [];
end

u_mask = repmat((masku./masku),1,1,40,size(vtran_mod,4));
v_mask = repmat((maskv./maskv),1,1,40,size(vtran_mod,4));

utran_mod = utran_mod .* u_mask;
vtran_mod = vtran_mod .* v_mask;


%% rain
for i = 1:length(near_p_1)
tran_v_1_pre(i,:)=sum(vtran_mod(near_p_1(i,1),near_p_1(i,2),:,:),[1 2 3],'omitnan');
end
for i = 1:length(near_p_2)
tran_u_2_pre(i,:)=sum(utran_mod(near_p_2(i,1),near_p_2(i,2),:,:),[1 2 3],'omitnan');
end
for i = 1:length(near_p_3)
tran_u_3_pre(i,:)=sum(utran_mod(near_p_3(i,1),near_p_3(i,2),:,:),[1 2 3],'omitnan');
end
% sec 4
for i = 1:length(p_4_v)
% tran_u_4(i,:)=sum(utran_mod(near_p_4_2(p_4_u(i),1),near_p_4_2(p_4_u(i),2),:,:),[1 2 3],'omitnan');    
tran_v_4_2_pre(i,:)=sum(vtran_mod(near_p_4_2(p_4_v(i),1),near_p_4_2(p_4_v(i),2),:,:),[1 2 3],'omitnan');    
end
for i = 1:length(p_4_u) 
tran_u_4_2_pre(i,:)=sum(utran_mod(near_p_4_2(p_4_u(i),1),near_p_4_2(p_4_u(i),2),:,:),[1 2 3],'omitnan'); 
end

for i = 1:length(near_p_4)   
tran_v_4_pre(i,:)=sum(vtran_mod(near_p_4(i,1),near_p_4(i,2),:,:),[1 2 3],'omitnan');    
end

% sec 5
for i = 1:length(p_5_v)  
tran_v_5_2_pre(i,:)=sum(vtran_mod(near_p_5_2(p_5_v(i),1),near_p_5_2(p_5_v(i),2)-p_5_v_1(i),:,:),[1 2 3],'omitnan');    
end
for i = 1:length(p_5_u) 
tran_u_5_2_pre(i,:)= sum(utran_mod(near_p_5_2(p_5_u(i),1),near_p_5_2(p_5_u(i),2),:,:),[1 2 3],'omitnan'); 
end

% sec 6
for i = 1:length(p_6_v)  
tran_v_6_2_pre(i,:)=sum(vtran_mod(near_p_6_2(p_6_v(i),1),near_p_6_2(p_6_v(i),2)-p_6_v_1(i),:,:),[1 2 3],'omitnan');    
end
for i = 1:length(p_6_u) 
tran_u_6_2_pre(i,:)= sum(utran_mod(near_p_6_2(p_6_u(i),1),near_p_6_2(p_6_u(i),2),:,:),[1 2 3],'omitnan'); 
end

if iy == 1982 
    tran_v_1 = tran_v_1_pre;
    tran_u_2 = tran_u_2_pre;

    tran_u_3 = tran_u_3_pre;
    tran_v_4_2 = tran_v_4_2_pre;
    tran_u_4_2 = tran_u_4_2_pre;
    tran_v_4 = tran_v_4_pre;

    tran_v_5_2 = tran_v_5_2_pre;
    tran_u_5_2 = tran_u_5_2_pre;
    tran_v_6_2 = tran_v_6_2_pre;
    tran_u_6_2 = tran_u_6_2_pre;

else
    tran_v_1 = cat(3,tran_v_1,tran_v_1_pre);
    tran_u_2 = cat(3,tran_u_2,tran_u_2_pre);

    tran_u_3 = cat(3,tran_u_3,tran_u_3_pre);
    tran_v_4_2 = cat(3,tran_v_4_2,tran_v_4_2_pre);
    tran_u_4_2 = cat(3,tran_u_4_2,tran_u_4_2_pre);
    tran_v_4 = cat(3,tran_v_4,tran_v_4_pre);

    tran_v_5_2 = cat(3,tran_v_5_2,tran_v_5_2_pre);
    tran_u_5_2 = cat(3,tran_u_5_2,tran_u_5_2_pre);
    tran_v_6_2 = cat(3,tran_v_6_2,tran_v_6_2_pre);
    tran_u_6_2 = cat(3,tran_u_6_2,tran_u_6_2_pre);
end

end

[raw txt]=xlsread('J:\fennel_NWP\Output\KoreaTsushima_Strait_transport_tidegauge(shin2022).xlsx','Sheet1','');

trans=raw(8:end,2:13)';
trans_yr=mean(trans,[1],'omitnan');
trans=reshape(trans,1,size(trans,1)*size(trans,2));
trans_1_m = reshape(squeeze(sum(tran_v_1,1)),1,12*39)./10^6;
trans_1_yr = mean(squeeze(sum(tran_v_1,1))./10^6,[1],'omitnan');

save('transport_1982to2020_noep_f.mat');


%% full time series
ontSizeTick = 13
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
plot(trans_1_m);
plot(trans,'r');
xlim([1 inf]);
xticks(1:24:length(1982:2020)*12); 
xticklabels([1982:2:2020]);
title({'KOREA Strait transport'});
ylabel('volume transport (Sv)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on'); 
xtickangle(45);
saveas(gcf,strcat('korea_strait_transport_1982to2020.png'),'png');
close; 


%% monthly
for i = 1:12
    ontSizeTick = 13
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
    plot(trans_1_m(i:12:end));
    plot(trans(i:12:end),'r');
xlim([1 inf]); ylim([1 4.5]);
xticks(1:2:length(1982:2020)); 
xticklabels([1982:2:2020]);
title(['KOREA Strait transport ',num2str(i),'month']);
ylabel('volume transport (Sv)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on'); 
xtickangle(45);
saveas(gcf,strcat('korea_strait_transport_1982to2020_',num2str(i),'-mon.png'),'png');
% close; 
end

%% yearly mean
 ontSizeTick = 13
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
plot(trans_1_yr);
plot(trans_yr,'r');
ylim([1 4.5]);
xticks(1:2:length(1982:2020)); 
xticklabels([1982:2:2020]);
title(['KOREA Strait transport yearly mean']);
ylabel('volume transport (Sv)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on'); 
xtickangle(45);
saveas(gcf,strcat('korea_strait_transport_1982to2020_yearly.png'),'png');

%% xy plot
% yearly mean
 ontSizeTick = 13
figPos = [0 0 5 4];
figure('Position',figPos*get(groot,'ScreenPixelsPerInch'),'paperunits','inches','paperposition',figPos);  hold on;
plot(trans,trans_1_m(1:length(trans)),'.');
p= polyfit(trans(1:441),trans_1_m(1:441),1);
y=polyval(p,0:0.5:5);
plot(0:0.5:5,y,'r');
r=corrcoef(trans(1:441),trans_1_m(1:441));
grid on; text(3.5,2.0,['r = ',num2str(r(1,2),'%0.2f')],'fontsize',15,'fontweight','bold');
xlim([1 4.5]); ylim([1 4.5]);
title(['KOREA Strait transport TG vs. model']);
ylabel('model volume transport (Sv)');
xlabel('TG volume transport (Sv)');
set(gca,'linewidth',1.5,'fontsize',fontSizeTick,'box','on'); grid(gca,'on'); 
saveas(gcf,strcat('korea_strait_transport_TG_vs_model_1982to2018.png'),'png');
