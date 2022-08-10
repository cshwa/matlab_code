%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 400_16_DO_surf. %%
%% first run the plot_KODC_4_40616_polynomial_partial
temp_s_nonan=find(isnan(temp_sur)==0);
yCalc = temp_sur(temp_s_nonan) .* y_s_do_1(1) + y_s_do_1(2);

temp_b_nonan=find(isnan(temp_bot)==0);
yCalc_b = temp_bot(temp_b_nonan) .* y_b_do_1(1) + y_b_do_1(2);

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= do_sur;
   temp_b= do_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2); 
plot(temp_s_nonan, yCalc,'*','color','r','linew',2);
%  plot(temp_b,'*','color','r','linew',2);
%  plot(temp_b_nonan, yCalc_b,'*','color','b','linew',2);
    
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
% temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

 plot(temp_s,'color','b');
%  plot(temp_b,'color','r');  
ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13) 
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);


% make input
% temp_s_nonan, yCalc
% temp_s_nonan(1:61)
clearvars temp_temp_s
temp_temp_s =do_sur;


temp_temp_s(temp_s_nonan(9:18))= yCalc(9:18);

temp_temp_s(isnan(temp_temp_s)) = interp1( t(~isnan(temp_temp_s)), temp_temp_s(~isnan(temp_temp_s)), t(isnan(temp_temp_s)) ); 
plot(temp_temp_s,'k')

do_re_s = temp_temp_s;
save('do_s_kodc_input_400_16.mat','tx_tick','do_re_s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 205_01_DO_surf. %%
%% first run the plot_KODC_4_20501_polynomial_partial
temp_s_nonan=find(isnan(temp_sur)==0);
yCalc = temp_sur(temp_s_nonan) .* y_s_do_1(1) + y_s_do_1(2);

temp_b_nonan=find(isnan(temp_bot)==0);
yCalc_b = temp_bot(temp_b_nonan) .* y_b_do_1(1) + y_b_do_1(2);

clearvars b1 yCalc1 X b yCalc2 nonan idx nanaxisx *_mth_f idx1 sig_txt *_mi2 *_bc
color_spec = ['k', 'r', 'b', 'g', 'c', 'm'];
gray = [128/255 128/255 128/255];
    clearvars temp_s temp_b
   temp_s= do_sur;
   temp_b= do_bot;
figure; hold on;
 plot(temp_s,'*','color','b','linew',2); 
plot(temp_s_nonan, yCalc,'*','color','r','linew',2);
%  plot(temp_b,'*','color','r','linew',2);
%  plot(temp_b_nonan, yCalc_b,'*','color','b','linew',2);
    
t=1:length(no3_sur);
%      temp_t_s(isnan(temp_t_s)) = interp1(t(~isnan(temp_t_s)), temp_t_s(~isnan(temp_t_s)), t(isnan(temp_t_s)) ); 
temp_s(isnan(temp_s)) = interp1( t(~isnan(temp_s)), temp_s(~isnan(temp_s)), t(isnan(temp_s)) ); 
% temp_b(isnan(temp_b)) = interp1( t(~isnan(temp_b)), temp_b(~isnan(temp_b)), t(isnan(temp_b)) ); 

interp_nh4_b = temp_b;
interp_nh4_s = temp_s;

 plot(temp_s,'color','b');
%  plot(temp_b,'color','r');  
ylim([-inf inf])
xlabel('time(year)','fontsize',13)
ylabel('do (mg/L)','fontsize',13) 
grid on; set(gca,'fontsize',13)
set(gca,'xtick',[tx_tick(1:1:11)]);
set(gca,'xlim',[tx_tick(1) tx_tick(11)]);
set(gca,'xticklabel',1980:1:1990);


% make input
% temp_s_nonan, yCalc
% temp_s_nonan(1:61)
clearvars temp_temp_s
temp_temp_s =do_sur;

temp_temp_s(temp_s_nonan(9:18))= yCalc(9:18);

temp_temp_s(isnan(temp_temp_s)) = interp1( t(~isnan(temp_temp_s)), temp_temp_s(~isnan(temp_temp_s)), t(isnan(temp_temp_s)) ); 
plot(temp_temp_s,'k')

do_re_s = temp_temp_s;
save('do_s_kodc_input_205_01.mat','tx_tick','do_re_s');


