close all; clc; clear;

% name_tag{sp_gy}
% ±¤¾çÇ×, ±¤¾ç4, ±¤¾ç3, ±¤¾ç2, ±¤¾ç1, ±¤¾ç5, ¿©¼ö2, ¿©¼ö3, ¿©¼ö1
% [res I]=sort([4,3,2,1,5]);
P_MW = 30.973762;
PO4_MW =94.971482;
N_MW = 14.006720;
NO3_MW = 62.005010;
NH4_MW = 18.038508;

yoonja=load('yoonjakangs_koem_data_monthly.mat');

% input the sigma
% sig = 3; %% sigma
sig = 2; %% sigma

% % station name tag
% % (for preventing the miss picking using tag from 1 col's name on excel) 
[res I]=sort([1,5,4,3,2,6,8,9,7]);
% sorting confirm
nn=['±¤¾çÇ×'; '±¤¾ç4'; '±¤¾ç3'; '±¤¾ç2'; '±¤¾ç1'; '±¤¾ç5'; '¿©¼ö2'; '¿©¼ö3'; '¿©¼ö1'];

% combining the tag and outter point excluding
name_tag = nn(I,:); 
size_tag = length(name_tag);

%% pick matched name with tag
% %port
% for i = 1:length(name_tag)
%    if  sum(strcmp(name_tag{i}, txt_matc_p)) ~= 0
%        indx{i} = find([strcmp(name_tag{i}, txt_matc_p)] == 1)     
%    end
% end


%% make 1997 to 2018 'yymm' form
k=0
for i = 1997:2019
    for j = 1:12
        k=k+1;
        ref_date{k,1} = [num2str(i) '-' num2str(j,'%02d')];
    end
end

%% make 1997 to 2018 'yymm' form
for j = 1:12
 ref_date_mm{j,1} = [num2str(j,'%02d')];
end

% spatial region index
do_sur_clim = yoonja.do_sur;
no3_sur_clim = yoonja.no3_sur;
temp_sur_clim = yoonja.temp_sur;
salt_sur_clim = yoonja.salt_sur;
chl_sur_clim = yoonja.chl_sur;
nh4_sur_clim = yoonja.nh4_sur;
        
po4_sur_clim = yoonja.po4_sur;
po4_bot_clim = yoonja.po4_bot;
        
do_bot_clim = yoonja.do_bot;
no3_bot_clim = yoonja.no3_bot;
temp_bot_clim = yoonja.temp_bot;
salt_bot_clim = yoonja.salt_bot;
chl_bot_clim = yoonja.chl_bot;
nh4_bot_clim = yoonja.nh4_bot;

%% pick only in the bay station (±¤¾çÇ×','±¤¾ç1','±¤¾ç2','±¤¾ç3');
do_sur_clim(1:5,:) = NaN;
no3_sur_clim(1:5,:) = NaN;
temp_sur_clim(1:5,:) = NaN;
salt_sur_clim(1:5,:) = NaN;
chl_sur_clim(1:5,:) = NaN;
nh4_sur_clim(1:5,:) = NaN;
        
po4_sur_clim(1:5,:) = NaN;
po4_bot_clim(1:5,:) = NaN;
        
do_bot_clim(1:5,:) = NaN;
no3_bot_clim(1:5,:) = NaN;
temp_bot_clim(1:5,:) = NaN;
salt_bot_clim(1:5,:) = NaN;
chl_bot_clim(1:5,:) = NaN;
nh4_bot_clim(1:5,:) = NaN;

sp_gy=size(yoonja.do_bot,1);
%% extract over 1~3sig for each st.
% for i = 1:length(name_tag)
    clearvars regime_*

    % sur
for i = 1:sp_gy
    clearvars idx1 tempo_data
    tempo_data= reshape(do_sur_clim(i,:),1,size(do_sur_clim(i,:),1)*size(do_sur_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_do=tempo_data;
    for sig_ind = 1:3
        do_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        do_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_do(find(tempo_data > do_up_bound(3)))=NaN;
    regime_do(find(tempo_data < do_down_bound(3)))=NaN;
    mtx_regime_do_3s(i,:) = reshape(regime_do,size(do_sur_clim(i,:),1),size(do_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(do_sur_clim(i,:),1,size(do_sur_clim(i,:),1)*size(do_sur_clim(i,:),2)); 
    regime_do=tempo_data;
    regime_do(find(tempo_data > do_up_bound(2)))=NaN;
    regime_do(find(tempo_data < do_down_bound(2)))=NaN;
    mtx_regime_do_2s(i,:) = reshape(regime_do,size(do_sur_clim(i,:),1),size(do_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(do_sur_clim(i,:),1,size(do_sur_clim(i,:),1)*size(do_sur_clim(i,:),2)); 
    regime_do=tempo_data;
    regime_do(find(tempo_data > do_up_bound(1)))=NaN;
    regime_do(find(tempo_data < do_down_bound(1)))=NaN;
    mtx_regime_do_1s(i,:) = reshape(regime_do,size(do_sur_clim(i,:),1),size(do_sur_clim(i,:),2));
    
    
   clearvars idx1 tempo_data
    tempo_data= reshape(no3_sur_clim(i,:),1,size(no3_sur_clim(i,:),1)*size(no3_sur_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_no3=tempo_data;
    for sig_ind = 1:3
        no3_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        no3_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_no3(find(tempo_data > no3_up_bound(3)))=NaN;
    regime_no3(find(tempo_data < no3_down_bound(3)))=NaN;
    mtx_regime_no3_3s(i,:) = reshape(regime_no3,size(no3_sur_clim(i,:),1),size(no3_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(no3_sur_clim(i,:),1,size(no3_sur_clim(i,:),1)*size(no3_sur_clim(i,:),2));
    regime_no3=tempo_data;
    regime_no3(find(tempo_data > no3_up_bound(2)))=NaN;
    regime_no3(find(tempo_data < no3_down_bound(2)))=NaN;
    mtx_regime_no3_2s(i,:) = reshape(regime_no3,size(no3_sur_clim(i,:),1),size(no3_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(no3_sur_clim(i,:),1,size(no3_sur_clim(i,:),1)*size(no3_sur_clim(i,:),2));
    regime_no3=tempo_data;
    regime_no3(find(tempo_data > no3_up_bound(1)))=NaN;
    regime_no3(find(tempo_data < no3_down_bound(1)))=NaN;
    mtx_regime_no3_1s(i,:) = reshape(regime_no3,size(no3_sur_clim(i,:),1),size(no3_sur_clim(i,:),2));
    

   clearvars idx1 tempo_data
    tempo_data= reshape(temp_sur_clim(i,:),1,size(temp_sur_clim(i,:),1)*size(temp_sur_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_temp=tempo_data;
    for sig_ind = 1:3
        temp_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        temp_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_temp(find(tempo_data > temp_up_bound(3)))=NaN;
    regime_temp(find(tempo_data < temp_down_bound(3)))=NaN;
    mtx_regime_temp_3s(i,:) = reshape(regime_temp,size(temp_sur_clim(i,:),1),size(temp_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(temp_sur_clim(i,:),1,size(temp_sur_clim(i,:),1)*size(temp_sur_clim(i,:),2));
    regime_temp=tempo_data;
    regime_temp(find(tempo_data > temp_up_bound(2)))=NaN;
    regime_temp(find(tempo_data < temp_down_bound(2)))=NaN;
    mtx_regime_temp_2s(i,:) = reshape(regime_temp,size(temp_sur_clim(i,:),1),size(temp_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(temp_sur_clim(i,:),1,size(temp_sur_clim(i,:),1)*size(temp_sur_clim(i,:),2));
    regime_temp=tempo_data;
    regime_temp(find(tempo_data > temp_up_bound(1)))=NaN;
    regime_temp(find(tempo_data < temp_down_bound(1)))=NaN;
    mtx_regime_temp_1s(i,:) = reshape(regime_temp,size(temp_sur_clim(i,:),1),size(temp_sur_clim(i,:),2));
    

    clearvars idx1 tempo_data
    tempo_data= reshape(salt_sur_clim(i,:),1,size(salt_sur_clim(i,:),1)*size(salt_sur_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_salt=tempo_data;
    for sig_ind = 1:3
        salt_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        salt_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_salt(find(tempo_data > salt_up_bound(3)))=NaN;
    regime_salt(find(tempo_data < salt_down_bound(3)))=NaN;
    mtx_regime_salt_3s(i,:) = reshape(regime_salt,size(salt_sur_clim(i,:),1),size(salt_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(salt_sur_clim(i,:),1,size(salt_sur_clim(i,:),1)*size(salt_sur_clim(i,:),2));
    regime_salt=tempo_data;
    regime_salt(find(tempo_data > salt_up_bound(2)))=NaN;
    regime_salt(find(tempo_data < salt_down_bound(2)))=NaN;
    mtx_regime_salt_2s(i,:) = reshape(regime_salt,size(salt_sur_clim(i,:),1),size(salt_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(salt_sur_clim(i,:),1,size(salt_sur_clim(i,:),1)*size(salt_sur_clim(i,:),2));
    regime_salt=tempo_data;
    regime_salt(find(tempo_data > salt_up_bound(1)))=NaN;
    regime_salt(find(tempo_data < salt_down_bound(1)))=NaN;
    mtx_regime_salt_1s(i,:) = reshape(regime_salt,size(salt_sur_clim(i,:),1),size(salt_sur_clim(i,:),2));
   


   clearvars idx1 tempo_data
    tempo_data= reshape(po4_sur_clim(i,:),1,size(po4_sur_clim(i,:),1)*size(po4_sur_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_po4=tempo_data;
    for sig_ind = 1:3
        po4_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        po4_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_po4(find(tempo_data > po4_up_bound(3)))=NaN;
    regime_po4(find(tempo_data < po4_down_bound(3)))=NaN;
    mtx_regime_po4_3s(i,:) = reshape(regime_po4,size(po4_sur_clim(i,:),1),size(po4_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(po4_sur_clim(i,:),1,size(po4_sur_clim(i,:),1)*size(po4_sur_clim(i,:),2));
    regime_po4=tempo_data;
    regime_po4(find(tempo_data > po4_up_bound(2)))=NaN;
    regime_po4(find(tempo_data < po4_down_bound(2)))=NaN;
    mtx_regime_po4_2s(i,:) = reshape(regime_po4,size(po4_sur_clim(i,:),1),size(po4_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(po4_sur_clim(i,:),1,size(po4_sur_clim(i,:),1)*size(po4_sur_clim(i,:),2));
    regime_po4=tempo_data;
    regime_po4(find(tempo_data > po4_up_bound(1)))=NaN;
    regime_po4(find(tempo_data < po4_down_bound(1)))=NaN;
    mtx_regime_po4_1s(i,:) = reshape(regime_po4,size(po4_sur_clim(i,:),1),size(po4_sur_clim(i,:),2));
    
    
%     figure; hold on;
%     for i = 1:sp_gy
%         plot(po4_sur_clim(i,:)./P_MW,'k.');
%     end
%     color_p={'r','g','b'};
%     for i = 1:3
%     yline(po4_up_bound(i)./P_MW,color_p{i}); yline(po4_down_bound(i)./P_MW,color_p{i});
%     end
%    ylim([0 inf])
%    
%    figure; hold on;
%     for i = 1:sp_gy
%         plot(mtx_regime_po4_1s(i,:)./P_MW,'k.');
%     end
%     color_p={'r','g','b'};
%     for i = 1:3
%     yline(po4_up_bound(i)./P_MW,color_p{i}); yline(po4_down_bound(i)./P_MW,color_p{i});
%     end
%    ylim([0 inf])
%    
%    mtx_1s = nanmean(mtx_regime_po4_1s,1);
%    mtx_2s = nanmean(mtx_regime_po4_2s,1);
%    mtx_3s = nanmean(mtx_regime_po4_3s,1);
%     figure; hold on;
%     plot(find(isnan(mtx_1s)==0),mtx_1s(~isnan(mtx_1s))./P_MW,'r');
%     plot(find(isnan(mtx_2s)==0),mtx_2s(~isnan(mtx_2s))./P_MW,'g');
%     plot(find(isnan(mtx_3s)==0),mtx_3s(~isnan(mtx_3s))./P_MW,'b');
%     
   
%     color_p={'r','g','b'};
%     for i = 1:3
%     yline(po4_up_bound(i)./P_MW,color_p{i}); yline(po4_down_bound(i)./P_MW,color_p{i});
%     end
%    ylim([0 inf])
    
    clearvars idx1 tempo_data
    tempo_data= reshape(nh4_sur_clim(i,:),1,size(nh4_sur_clim(i,:),1)*size(nh4_sur_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_nh4=tempo_data;
    for sig_ind = 1:3
        nh4_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        nh4_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_nh4(find(tempo_data > nh4_up_bound(3)))=NaN;
    regime_nh4(find(tempo_data < nh4_down_bound(3)))=NaN;
    mtx_regime_nh4_3s(i,:) = reshape(regime_nh4,size(nh4_sur_clim(i,:),1),size(nh4_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(nh4_sur_clim(i,:),1,size(nh4_sur_clim(i,:),1)*size(nh4_sur_clim(i,:),2));
    regime_nh4=tempo_data;
    regime_nh4(find(tempo_data > nh4_up_bound(2)))=NaN;
    regime_nh4(find(tempo_data < nh4_down_bound(2)))=NaN;
    mtx_regime_nh4_2s(i,:) = reshape(regime_nh4,size(nh4_sur_clim(i,:),1),size(nh4_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(nh4_sur_clim(i,:),1,size(nh4_sur_clim(i,:),1)*size(nh4_sur_clim(i,:),2));
    regime_nh4=tempo_data;
    regime_nh4(find(tempo_data > nh4_up_bound(1)))=NaN;
    regime_nh4(find(tempo_data < nh4_down_bound(1)))=NaN;
    mtx_regime_nh4_1s(i,:) = reshape(regime_nh4,size(nh4_sur_clim(i,:),1),size(nh4_sur_clim(i,:),2));
    
    
   clearvars idx1 tempo_data
    tempo_data= reshape(chl_sur_clim(i,:),1,size(chl_sur_clim(i,:),1)*size(chl_sur_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_chl=tempo_data;
    for sig_ind = 1:3
        chl_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        chl_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_chl(find(tempo_data > chl_up_bound(3)))=NaN;
    regime_chl(find(tempo_data < chl_down_bound(3)))=NaN;
    mtx_regime_chl_3s(i,:) = reshape(regime_chl,size(chl_sur_clim(i,:),1),size(chl_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(chl_sur_clim(i,:),1,size(chl_sur_clim(i,:),1)*size(chl_sur_clim(i,:),2));
    regime_chl=tempo_data;
    regime_chl(find(tempo_data > chl_up_bound(2)))=NaN;
    regime_chl(find(tempo_data < chl_down_bound(2)))=NaN;
    mtx_regime_chl_2s(i,:) = reshape(regime_chl,size(chl_sur_clim(i,:),1),size(chl_sur_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(chl_sur_clim(i,:),1,size(chl_sur_clim(i,:),1)*size(chl_sur_clim(i,:),2));
    regime_chl=tempo_data;
    regime_chl(find(tempo_data > chl_up_bound(1)))=NaN;
    regime_chl(find(tempo_data < chl_down_bound(1)))=NaN;
    mtx_regime_chl_1s(i,:) = reshape(regime_chl,size(chl_sur_clim(i,:),1),size(chl_sur_clim(i,:),2));
    
     %% bot
    clearvars idx1 tempo_data
    tempo_data= reshape(do_bot_clim(i,:),1,size(do_bot_clim(i,:),1)*size(do_bot_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_do_b=tempo_data;
    for sig_ind = 1:3
        do_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        do_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_do_b(find(tempo_data > do_b_up_bound(3)))=NaN;
    regime_do_b(find(tempo_data < do_b_down_bound(3)))=NaN;
    mtx_regime_do_b_3s(i,:) = reshape(regime_do_b,size(do_bot_clim(i,:),1),size(do_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(do_bot_clim(i,:),1,size(do_bot_clim(i,:),1)*size(do_bot_clim(i,:),2));
    regime_do_b=tempo_data;
    regime_do_b(find(tempo_data > do_b_up_bound(2)))=NaN;
    regime_do_b(find(tempo_data < do_b_down_bound(2)))=NaN;
    mtx_regime_do_b_2s(i,:) = reshape(regime_do_b,size(do_bot_clim(i,:),1),size(do_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(do_bot_clim(i,:),1,size(do_bot_clim(i,:),1)*size(do_bot_clim(i,:),2));
    regime_do_b=tempo_data;
    regime_do_b(find(tempo_data > do_b_up_bound(1)))=NaN;
    regime_do_b(find(tempo_data < do_b_down_bound(1)))=NaN;
    mtx_regime_do_b_1s(i,:) = reshape(regime_do_b,size(do_bot_clim(i,:),1),size(do_bot_clim(i,:),2));
        
    
  clearvars idx1 tempo_data
    tempo_data= reshape(no3_bot_clim(i,:),1,size(no3_bot_clim(i,:),1)*size(no3_bot_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_no3_b=tempo_data;
    for sig_ind = 1:3
        no3_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        no3_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_no3_b(find(tempo_data > no3_b_up_bound(3)))=NaN;
    regime_no3_b(find(tempo_data < no3_b_down_bound(3)))=NaN;
    mtx_regime_no3_b_3s(i,:) = reshape(regime_no3_b,size(no3_bot_clim(i,:),1),size(no3_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(no3_bot_clim(i,:),1,size(no3_bot_clim(i,:),1)*size(no3_bot_clim(i,:),2));
    regime_no3_b=tempo_data;
    regime_no3_b(find(tempo_data > no3_b_up_bound(2)))=NaN;
    regime_no3_b(find(tempo_data < no3_b_down_bound(2)))=NaN;
    mtx_regime_no3_b_2s(i,:) = reshape(regime_no3_b,size(no3_bot_clim(i,:),1),size(no3_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(no3_bot_clim(i,:),1,size(no3_bot_clim(i,:),1)*size(no3_bot_clim(i,:),2));
    regime_no3_b=tempo_data;
    regime_no3_b(find(tempo_data > no3_b_up_bound(1)))=NaN;
    regime_no3_b(find(tempo_data < no3_b_down_bound(1)))=NaN;
    mtx_regime_no3_b_1s(i,:) = reshape(regime_no3_b,size(no3_bot_clim(i,:),1),size(no3_bot_clim(i,:),2));
    

   clearvars idx1 tempo_data
    tempo_data= reshape(temp_bot_clim(i,:),1,size(temp_bot_clim(i,:),1)*size(temp_bot_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_temp_b=tempo_data;
    for sig_ind = 1:3
        temp_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        temp_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_temp_b(find(tempo_data > temp_b_up_bound(3)))=NaN;
    regime_temp_b(find(tempo_data < temp_b_down_bound(3)))=NaN;
    mtx_regime_temp_b_3s(i,:) = reshape(regime_temp_b,size(temp_bot_clim(i,:),1),size(temp_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(temp_bot_clim(i,:),1,size(temp_bot_clim(i,:),1)*size(temp_bot_clim(i,:),2));
    regime_temp_b=tempo_data;
    regime_temp_b(find(tempo_data > temp_b_up_bound(2)))=NaN;
    regime_temp_b(find(tempo_data < temp_b_down_bound(2)))=NaN;
    mtx_regime_temp_b_2s(i,:) = reshape(regime_temp_b,size(temp_bot_clim(i,:),1),size(temp_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(temp_bot_clim(i,:),1,size(temp_bot_clim(i,:),1)*size(temp_bot_clim(i,:),2));
    regime_temp_b=tempo_data;
    regime_temp_b(find(tempo_data > temp_b_up_bound(1)))=NaN;
    regime_temp_b(find(tempo_data < temp_b_down_bound(1)))=NaN;
    mtx_regime_temp_b_1s(i,:) = reshape(regime_temp_b,size(temp_bot_clim(i,:),1),size(temp_bot_clim(i,:),2));
    

   clearvars idx1 tempo_data
    tempo_data= reshape(salt_bot_clim(i,:),1,size(salt_bot_clim(i,:),1)*size(salt_bot_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_salt_b=tempo_data;
    for sig_ind = 1:3
        salt_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        salt_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_salt_b(find(tempo_data > salt_b_up_bound(3)))=NaN;
    regime_salt_b(find(tempo_data < salt_b_down_bound(3)))=NaN;
    mtx_regime_salt_b_3s(i,:) = reshape(regime_salt_b,size(salt_bot_clim(i,:),1),size(salt_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(salt_bot_clim(i,:),1,size(salt_bot_clim(i,:),1)*size(salt_bot_clim(i,:),2));
    regime_salt_b=tempo_data;
    regime_salt_b(find(tempo_data > salt_b_up_bound(2)))=NaN;
    regime_salt_b(find(tempo_data < salt_b_down_bound(2)))=NaN;
    mtx_regime_salt_b_2s(i,:) = reshape(regime_salt_b,size(salt_bot_clim(i,:),1),size(salt_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(salt_bot_clim(i,:),1,size(salt_bot_clim(i,:),1)*size(salt_bot_clim(i,:),2));
    regime_salt_b=tempo_data;
    regime_salt_b(find(tempo_data > salt_b_up_bound(1)))=NaN;
    regime_salt_b(find(tempo_data < salt_b_down_bound(1)))=NaN;
    mtx_regime_salt_b_1s(i,:) = reshape(regime_salt_b,size(salt_bot_clim(i,:),1),size(salt_bot_clim(i,:),2));
    

    clearvars idx1 tempo_data
    tempo_data= reshape(po4_bot_clim,1,size(po4_bot_clim,1)*size(po4_bot_clim,2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_po4_b=tempo_data;
    for sig_ind = 1:3
        po4_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        po4_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_po4_b(find(tempo_data > po4_b_up_bound(3)))=NaN;
    regime_po4_b(find(tempo_data < po4_b_down_bound(3)))=NaN;
    mtx_regime_po4_b_3s = reshape(regime_po4_b,size(po4_bot_clim,1),size(po4_bot_clim,2));
    
    clearvars idx1 tempo_data
    tempo_data= reshape(po4_bot_clim(i,:),1,size(po4_bot_clim(i,:),1)*size(po4_bot_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_po4_b=tempo_data;
    for sig_ind = 1:3
        po4_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        po4_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_po4_b(find(tempo_data > po4_b_up_bound(3)))=NaN;
    regime_po4_b(find(tempo_data < po4_b_down_bound(3)))=NaN;
    mtx_regime_po4_b_3s(i,:) = reshape(regime_po4_b,size(po4_bot_clim(i,:),1),size(po4_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(po4_bot_clim(i,:),1,size(po4_bot_clim(i,:),1)*size(po4_bot_clim(i,:),2));
    regime_po4_b=tempo_data;
    regime_po4_b(find(tempo_data > po4_b_up_bound(2)))=NaN;
    regime_po4_b(find(tempo_data < po4_b_down_bound(2)))=NaN;
    mtx_regime_po4_b_2s(i,:) = reshape(regime_po4_b,size(po4_bot_clim(i,:),1),size(po4_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(po4_bot_clim(i,:),1,size(po4_bot_clim(i,:),1)*size(po4_bot_clim(i,:),2));
    regime_po4_b=tempo_data;
    regime_po4_b(find(tempo_data > po4_b_up_bound(1)))=NaN;
    regime_po4_b(find(tempo_data < po4_b_down_bound(1)))=NaN;
    mtx_regime_po4_b_1s(i,:) = reshape(regime_po4_b,size(po4_bot_clim(i,:),1),size(po4_bot_clim(i,:),2));
    

    clearvars idx1 tempo_data
    tempo_data= reshape(nh4_bot_clim(i,:),1,size(nh4_bot_clim(i,:),1)*size(nh4_bot_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_nh4_b=tempo_data;
    for sig_ind = 1:3
        nh4_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        nh4_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_nh4_b(find(tempo_data > nh4_b_up_bound(3)))=NaN;
    regime_nh4_b(find(tempo_data < nh4_b_down_bound(3)))=NaN;
    mtx_regime_nh4_b_3s(i,:) = reshape(regime_nh4_b,size(nh4_bot_clim(i,:),1),size(nh4_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(nh4_bot_clim(i,:),1,size(nh4_bot_clim(i,:),1)*size(nh4_bot_clim(i,:),2));
    regime_nh4_b=tempo_data;
    regime_nh4_b(find(tempo_data > nh4_b_up_bound(2)))=NaN;
    regime_nh4_b(find(tempo_data < nh4_b_down_bound(2)))=NaN;
    mtx_regime_nh4_b_2s(i,:) = reshape(regime_nh4_b,size(nh4_bot_clim(i,:),1),size(nh4_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(nh4_bot_clim(i,:),1,size(nh4_bot_clim(i,:),1)*size(nh4_bot_clim(i,:),2));
    regime_nh4_b=tempo_data;
    regime_nh4_b(find(tempo_data > nh4_b_up_bound(1)))=NaN;
    regime_nh4_b(find(tempo_data < nh4_b_down_bound(1)))=NaN;
    mtx_regime_nh4_b_1s(i,:) = reshape(regime_nh4_b,size(nh4_bot_clim(i,:),1),size(nh4_bot_clim(i,:),2));
    
    
   clearvars idx1 tempo_data
    tempo_data= reshape(chl_bot_clim(i,:),1,size(chl_bot_clim(i,:),1)*size(chl_bot_clim(i,:),2)); 
    idx1 = find(isnan(tempo_data) == 0);
    regime_chl_b=tempo_data;
    for sig_ind = 1:3
        chl_b_up_bound(sig_ind)= mean(tempo_data(idx1)) + sig_ind*std(tempo_data(idx1));
        chl_b_down_bound(sig_ind)= mean(tempo_data(idx1)) - sig_ind*std(tempo_data(idx1));
    end
    regime_chl_b(find(tempo_data > chl_b_up_bound(3)))=NaN;
    regime_chl_b(find(tempo_data < chl_b_down_bound(3)))=NaN;
    mtx_regime_chl_b_3s(i,:) = reshape(regime_chl_b,size(chl_bot_clim(i,:),1),size(chl_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(chl_bot_clim(i,:),1,size(chl_bot_clim(i,:),1)*size(chl_bot_clim(i,:),2));
    regime_chl_b=tempo_data;
    regime_chl_b(find(tempo_data > chl_b_up_bound(2)))=NaN;
    regime_chl_b(find(tempo_data < chl_b_down_bound(2)))=NaN;
    mtx_regime_chl_b_2s(i,:) = reshape(regime_chl_b,size(chl_bot_clim(i,:),1),size(chl_bot_clim(i,:),2));
    
    clearvars tempo_data regime_*
    tempo_data= reshape(chl_bot_clim(i,:),1,size(chl_bot_clim(i,:),1)*size(chl_bot_clim(i,:),2));
    regime_chl_b=tempo_data;
    regime_chl_b(find(tempo_data > chl_b_up_bound(1)))=NaN;
    regime_chl_b(find(tempo_data < chl_b_down_bound(1)))=NaN;
    mtx_regime_chl_b_1s(i,:) = reshape(regime_chl_b,size(chl_bot_clim(i,:),1),size(chl_bot_clim(i,:),2));
    
   
end
    
%% extract regime mean
% regime_do = regime_do - regm_do;
% regime_no3 = regime_no3 - regm_no3;
% regime_temp = regime_temp - regm_temp;
% regime_salt = regime_salt - regm_salt;
% 
% regime_do_b = regime_do_b - regm_do_b;
% regime_no3_b = regime_no3_b - regm_no3_b;
% regime_temp_b = regime_temp_b - regm_temp_b;
% regime_salt_b = regime_salt_b - regm_salt_b;
% end

clearvars regime_*
regime_do_1s = mtx_regime_do_1s;
regime_no3_1s = mtx_regime_no3_1s;
regime_temp_1s = mtx_regime_temp_1s;
regime_salt_1s = mtx_regime_salt_1s;
regime_nh4_1s = mtx_regime_nh4_1s;
regime_chl_1s = mtx_regime_chl_1s;

regime_do_b_1s = mtx_regime_do_b_1s;
regime_no3_b_1s = mtx_regime_no3_b_1s;
regime_temp_b_1s = mtx_regime_temp_b_1s;
regime_salt_b_1s = mtx_regime_salt_b_1s;
regime_nh4_b_1s = mtx_regime_nh4_b_1s;
regime_chl_b_1s = mtx_regime_chl_b_1s;

regime_po4_1s = mtx_regime_po4_1s;
regime_po4_b_1s = mtx_regime_po4_b_1s;

regime_do_2s = mtx_regime_do_2s;
regime_no3_2s = mtx_regime_no3_2s;
regime_temp_2s = mtx_regime_temp_2s;
regime_salt_2s = mtx_regime_salt_2s;
regime_nh4_2s = mtx_regime_nh4_2s;
regime_chl_2s = mtx_regime_chl_2s;

regime_do_b_2s = mtx_regime_do_b_2s;
regime_no3_b_2s = mtx_regime_no3_b_2s;
regime_temp_b_2s = mtx_regime_temp_b_2s;
regime_salt_b_2s = mtx_regime_salt_b_2s;
regime_nh4_b_2s = mtx_regime_nh4_b_2s;
regime_chl_b_2s = mtx_regime_chl_b_2s;

regime_po4_2s = mtx_regime_po4_2s;
regime_po4_b_2s = mtx_regime_po4_b_2s;

regime_do_3s = mtx_regime_do_3s;
regime_no3_3s = mtx_regime_no3_3s;
regime_temp_3s = mtx_regime_temp_3s;
regime_salt_3s = mtx_regime_salt_3s;
regime_nh4_3s = mtx_regime_nh4_3s;
regime_chl_3s = mtx_regime_chl_3s;

regime_do_b_3s = mtx_regime_do_b_3s;
regime_no3_b_3s = mtx_regime_no3_b_3s;
regime_temp_b_3s = mtx_regime_temp_b_3s;
regime_salt_b_3s = mtx_regime_salt_b_3s;
regime_nh4_b_3s = mtx_regime_nh4_b_3s;
regime_chl_b_3s = mtx_regime_chl_b_3s;

regime_po4_3s = mtx_regime_po4_3s;
regime_po4_b_3s = mtx_regime_po4_b_3s;

% save('koem_timeseires_monthly_gy_only_2sig_v2.mat');
save(['koem_timeseires_monthly_gy_only_1to3sig_v2_each_st_out_bay.mat']);
return



% right_here
%% make linear coeff. and save file

%% salt
color_pic = lines(size(regime_salt_yr,1));
marker_sty = {'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<',...
    'o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<','o','+','x','^','>','h','p','s','d','.','*','v','<'};
xp = 1:22;
j=0
figure; hold on;
for i = 1:size(regime_salt_yr,1)
  clearvars reg_data_salt xp_w_salt pf_w_salt
reg_data_salt = regime_salt_yr(i,:);
if isnan(reg_data_salt(1)) == 0 
    j = j+1;
    xp_w_salt = find(isnan(reg_data_salt)==0);
    pf_w_salt = polyfit(xp_w_salt,reg_data_salt(xp_w_salt),1);
    yp_w_salt(i,:) = polyval(pf_w_salt,xp);
    scatter(1:22,regime_salt_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_salt(i,:),'color',color_pic(i,:));
    coeff_salt(i,:) = pf_w_salt;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('salt (mg/m^3)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM-Ç¥Ãþ¿°ºÐ ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])


%% no3
j=0
figure; hold on;
for i = 1:size(regime_no3_yr,1)
  clearvars reg_data_no3 xp_w_no3 pf_w_no3
reg_data_no3 = regime_no3_yr(i,:);
if isnan(reg_data_no3(1)) == 0 
    j = j+1;
    xp_w_no3 = find(isnan(reg_data_no3)==0);
    pf_w_no3 = polyfit(xp_w_no3,reg_data_no3(xp_w_no3),1);
    yp_w_no3(i,:) = polyval(pf_w_no3,xp);
    scatter(1:22,regime_no3_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_no3(i,:),'color',color_pic(i,:));
    coeff_no3(i,:) = pf_w_no3;
    temp_case(j) = i;
end
hold on
end
xlabel('time(year)','fontsize',13)
ylabel('NO3 (ug/L)','fontsize',13)
set(gca,'xtick',[1:2:22]);
set(gca,'xlim',[1 22]);
set(gca,'xticklabel',1997:2:2018);
title('KOEM°üÃø-Ç¥ÃþÁú»ê¿° ¿¬Æò±Õ','fontsize',13)
grid on
set(gca,'fontsize',13)
% ylim([32 35])
% legend('205-01','205-02','205-03','205-04','205-05')

%temp
j=0
figure; hold on;
for i = 1:size(regime_temp_yr,1)
  clearvars reg_data_temp xp_w_temp pf_w_temp
reg_data_temp = regime_temp_yr(i,:);
if isnan(reg_data_temp(1)) == 0 
    j = j+1;
    xp_w_temp = find(isnan(reg_data_temp)==0);
    pf_w_temp = polyfit(xp_w_temp,reg_data_temp(xp_w_temp),1);
    yp_w_temp(i,:) = polyval(pf_w_temp,xp);
    scatter(1:22,regime_temp_yr(i,:),marker_sty{i},'MarkerEdgeColor',color_pic(i,:));
    plot(1:22, yp_w_temp(i,:),'color',color_pic(i,:));
    coeff_temp(i,:) = pf_w_temp;
    temp_case(j) = i;
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

save('kodc_just_mean_yearly_linear_trend.mat','-v7.3');

filename=['KOEM_chl_bot_' num2str(i) 'mth']; 

print('-dpng',filename); 


close all

    % save('KOEM_name_tag.mat','name_tag');