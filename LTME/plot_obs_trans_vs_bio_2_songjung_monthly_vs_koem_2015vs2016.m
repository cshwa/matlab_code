
%% Kd vs. CHL

% year -1 11mth to 12mth
k=20
% surf
d_t_2008=nanmean(koem_temp_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_t_2008=nanstd(koem_temp_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_s_2008=nanmean(koem_salt_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_s_2008=nanstd(koem_salt_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_nh4_2008=nanmean(koem_nh4_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_nh4_2008=nanstd(koem_nh4_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_no3_2008=nanmean(koem_no3_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_no3_2008=nanstd(koem_no3_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_po4_2008=nanmean(koem_po4_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_po4_2008=nanstd(koem_po4_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_chl_2008=nanmean(koem_chl_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_chl_2008=nanstd(koem_chl_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_din_2008=nanmean(koem_din_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_din_2008=nanstd(koem_din_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_ss_2008=nanmean(koem_ss_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_ss_2008=nanstd(koem_ss_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_secchi_2008=nanmean(koem_secchi_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_secchi_2008=nanstd(koem_secchi_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_si_2008=nanmean(koem_si_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_si_2008=nanstd(koem_si_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_do_2008=nanmean(koem_do_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_do_2008=nanstd(koem_do_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
%bot
d_b_t_2008=nanmean(koem_temp_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_t_2008=nanstd(koem_temp_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_s_2008=nanmean(koem_salt_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_s_2008=nanstd(koem_salt_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_nh4_2008=nanmean(koem_nh4_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_nh4_2008=nanstd(koem_nh4_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_no3_2008=nanmean(koem_no3_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_no3_2008=nanstd(koem_no3_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_po4_2008=nanmean(koem_po4_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_po4_2008=nanstd(koem_po4_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_chl_2008=nanmean(koem_chl_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_chl_2008=nanstd(koem_chl_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_din_2008=nanmean(koem_din_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_din_2008=nanstd(koem_din_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_ss_2008=nanmean(koem_ss_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_ss_2008=nanstd(koem_ss_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_si_2008=nanmean(koem_si_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_si_2008=nanstd(koem_si_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
d_b_do_2008=nanmean(koem_do_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);
std_b_do_2008=nanstd(koem_do_b_raw(sp_9p_st,(k-1)*12+1 -2:12*k),1);

l=25
% surf
d_t_2013=nanmean(koem_temp_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_t_2013=nanstd(koem_temp_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_s_2013=nanmean(koem_salt_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_s_2013=nanstd(koem_salt_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_nh4_2013=nanmean(koem_nh4_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_nh4_2013=nanstd(koem_nh4_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_no3_2013=nanmean(koem_no3_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_no3_2013=nanstd(koem_no3_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_po4_2013=nanmean(koem_po4_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_po4_2013=nanstd(koem_po4_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_chl_2013=nanmean(koem_chl_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_chl_2013=nanstd(koem_chl_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_din_2013=nanmean(koem_din_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_din_2013=nanstd(koem_din_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_ss_2013=nanmean(koem_ss_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_ss_2013=nanstd(koem_ss_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_secchi_2013=nanmean(koem_secchi_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_secchi_2013=nanstd(koem_secchi_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_si_2013=nanmean(koem_si_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_si_2013=nanstd(koem_si_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_do_2013=nanmean(koem_do_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_do_2013=nanstd(koem_do_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
%bot
d_b_t_2013=nanmean(koem_temp_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_t_2013=nanstd(koem_temp_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_s_2013=nanmean(koem_salt_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_s_2013=nanstd(koem_salt_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_nh4_2013=nanmean(koem_nh4_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_nh4_2013=nanstd(koem_nh4_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_no3_2013=nanmean(koem_no3_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_no3_2013=nanstd(koem_no3_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_po4_2013=nanmean(koem_po4_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_po4_2013=nanstd(koem_po4_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_chl_2013=nanmean(koem_chl_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_chl_2013=nanstd(koem_chl_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_din_2013=nanmean(koem_din_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_din_2013=nanstd(koem_din_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_ss_2013=nanmean(koem_ss_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_ss_2013=nanstd(koem_ss_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_si_2013=nanmean(koem_si_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_si_2013=nanstd(koem_si_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
d_b_do_2013=nanmean(koem_do_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);
std_b_do_2013=nanstd(koem_do_b_raw(sp_9p_st,(l-1)*12+1 -2:12*l),1);

i=27
% surf
d_t_2015=nanmean(koem_temp_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_t_2015=nanstd(koem_temp_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_s_2015=nanmean(koem_salt_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_s_2015=nanstd(koem_salt_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_nh4_2015=nanmean(koem_nh4_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_nh4_2015=nanstd(koem_nh4_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_no3_2015=nanmean(koem_no3_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_no3_2015=nanstd(koem_no3_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_po4_2015=nanmean(koem_po4_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_po4_2015=nanstd(koem_po4_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_chl_2015=nanmean(koem_chl_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_chl_2015=nanstd(koem_chl_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_din_2015=nanmean(koem_din_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_din_2015=nanstd(koem_din_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_ss_2015=nanmean(koem_ss_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_ss_2015=nanstd(koem_ss_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_secchi_2015=nanmean(koem_secchi_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_secchi_2015=nanstd(koem_secchi_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_si_2015=nanmean(koem_si_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_si_2015=nanstd(koem_si_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_do_2015=nanmean(koem_do_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_do_2015=nanstd(koem_do_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
%bot
d_b_t_2015=nanmean(koem_temp_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_t_2015=nanstd(koem_temp_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_s_2015=nanmean(koem_salt_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_s_2015=nanstd(koem_salt_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_nh4_2015=nanmean(koem_nh4_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_nh4_2015=nanstd(koem_nh4_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_no3_2015=nanmean(koem_no3_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_no3_2015=nanstd(koem_no3_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_po4_2015=nanmean(koem_po4_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_po4_2015=nanstd(koem_po4_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_chl_2015=nanmean(koem_chl_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_chl_2015=nanstd(koem_chl_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_din_2015=nanmean(koem_din_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_din_2015=nanstd(koem_din_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_ss_2015=nanmean(koem_ss_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_ss_2015=nanstd(koem_ss_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_si_2015=nanmean(koem_si_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_si_2015=nanstd(koem_si_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
d_b_do_2015=nanmean(koem_do_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);
std_b_do_2015=nanstd(koem_do_b_raw(sp_9p_st,(i-1)*12+1 -2:12*i),1);

j=28
%surf
d_t_2016=nanmean(koem_temp_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_t_2016=nanstd(koem_temp_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_s_2016=nanmean(koem_salt_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_s_2016=nanstd(koem_salt_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_nh4_2016=nanmean(koem_nh4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_nh4_2016=nanstd(koem_nh4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_no3_2016=nanmean(koem_no3_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_no3_2016=nanstd(koem_no3_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_po4_2016=nanmean(koem_po4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_po4_2016=nanstd(koem_po4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_chl_2016=nanmean(koem_chl_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_chl_2016=nanstd(koem_chl_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_din_2016=nanmean(koem_din_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_din_2016=nanstd(koem_din_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_ss_2016=nanmean(koem_ss_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_ss_2016=nanstd(koem_ss_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_secchi_2016=nanmean(koem_secchi_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_secchi_2016=nanstd(koem_secchi_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_si_2016=nanmean(koem_si_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_si_2016=nanstd(koem_si_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_do_2016=nanmean(koem_do_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_do_2016=nanstd(koem_do_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
%bot
d_b_t_2016=nanmean(koem_temp_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_t_2016=nanstd(koem_temp_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_s_2016=nanmean(koem_salt_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_s_2016=nanstd(koem_salt_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_nh4_2016=nanmean(koem_nh4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_nh4_2016=nanstd(koem_nh4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_no3_2016=nanmean(koem_no3_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_no3_2016=nanstd(koem_no3_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_po4_2016=nanmean(koem_po4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_po4_2016=nanstd(koem_po4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_chl_2016=nanmean(koem_chl_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_chl_2016=nanstd(koem_chl_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_din_2016=nanmean(koem_din_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_din_2016=nanstd(koem_din_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_ss_2016=nanmean(koem_ss_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_ss_2016=nanstd(koem_ss_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_si_2016=nanmean(koem_si_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_si_2016=nanstd(koem_si_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_do_2016=nanmean(koem_do_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_do_2016=nanstd(koem_do_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);


j=29
%surf
d_t_2017=nanmean(koem_temp_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_t_2017=nanstd(koem_temp_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_s_2017=nanmean(koem_salt_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_s_2017=nanstd(koem_salt_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_nh4_2017=nanmean(koem_nh4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_nh4_2017=nanstd(koem_nh4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_no3_2017=nanmean(koem_no3_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_no3_2017=nanstd(koem_no3_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_po4_2017=nanmean(koem_po4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_po4_2017=nanstd(koem_po4_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_chl_2017=nanmean(koem_chl_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_chl_2017=nanstd(koem_chl_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_din_2017=nanmean(koem_din_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_din_2017=nanstd(koem_din_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_ss_2017=nanmean(koem_ss_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_ss_2017=nanstd(koem_ss_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_secchi_2017=nanmean(koem_secchi_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_secchi_2017=nanstd(koem_secchi_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_si_2017=nanmean(koem_si_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_si_2017=nanstd(koem_si_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_do_2017=nanmean(koem_do_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_do_2017=nanstd(koem_do_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
%bot
d_b_t_2017=nanmean(koem_temp_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_t_2017=nanstd(koem_temp_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_s_2017=nanmean(koem_salt_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_s_2017=nanstd(koem_salt_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_nh4_2017=nanmean(koem_nh4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_nh4_2017=nanstd(koem_nh4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_no3_2017=nanmean(koem_no3_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_no3_2017=nanstd(koem_no3_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_po4_2017=nanmean(koem_po4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_po4_2017=nanstd(koem_po4_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_chl_2017=nanmean(koem_chl_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_chl_2017=nanstd(koem_chl_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_din_2017=nanmean(koem_din_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_din_2017=nanstd(koem_din_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_ss_2017=nanmean(koem_ss_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_ss_2017=nanstd(koem_ss_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_si_2017=nanmean(koem_si_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_si_2017=nanstd(koem_si_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
d_b_do_2017=nanmean(koem_do_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);
std_b_do_2017=nanstd(koem_do_b_raw(sp_9p_st,(j-1)*12+1 -2:12*j),1);





% i=27
% % surf
% d_t_2015=nanmean(koem_temp_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_t_2015=nanstd(koem_temp_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_s_2015=nanmean(koem_salt_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_s_2015=nanstd(koem_salt_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_nh4_2015=nanmean(koem_nh4_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_nh4_2015=nanstd(koem_nh4_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_no3_2015=nanmean(koem_no3_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_no3_2015=nanstd(koem_no3_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_po4_2015=nanmean(koem_po4_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_po4_2015=nanstd(koem_po4_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_chl_2015=nanmean(koem_chl_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_chl_2015=nanstd(koem_chl_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_din_2015=nanmean(koem_din_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_din_2015=nanstd(koem_din_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_ss_2015=nanmean(koem_ss_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_ss_2015=nanstd(koem_ss_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_secchi_2015=nanmean(koem_secchi_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_secchi_2015=nanstd(koem_secchi_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_si_2015=nanmean(koem_si_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_si_2015=nanstd(koem_si_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_do_2015=nanmean(koem_do_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_do_2015=nanstd(koem_do_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% %bot
% d_b_t_2015=nanmean(koem_temp_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_t_2015=nanstd(koem_temp_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_s_2015=nanmean(koem_salt_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_s_2015=nanstd(koem_salt_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_nh4_2015=nanmean(koem_nh4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_nh4_2015=nanstd(koem_nh4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_no3_2015=nanmean(koem_no3_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_no3_2015=nanstd(koem_no3_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_po4_2015=nanmean(koem_po4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_po4_2015=nanstd(koem_po4_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_chl_2015=nanmean(koem_chl_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_chl_2015=nanstd(koem_chl_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_din_2015=nanmean(koem_din_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_din_2015=nanstd(koem_din_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_ss_2015=nanmean(koem_ss_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_ss_2015=nanstd(koem_ss_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_si_2015=nanmean(koem_si_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_si_2015=nanstd(koem_si_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% d_b_do_2015=nanmean(koem_do_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% std_b_do_2015=nanstd(koem_do_b_raw(sp_9p_st,(i-1)*12+1:12*i),1);
% 
% j=28
% %surf
% d_t_2016=nanmean(koem_temp_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_t_2016=nanstd(koem_temp_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_s_2016=nanmean(koem_salt_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_s_2016=nanstd(koem_salt_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_nh4_2016=nanmean(koem_nh4_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_nh4_2016=nanstd(koem_nh4_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_no3_2016=nanmean(koem_no3_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_no3_2016=nanstd(koem_no3_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_po4_2016=nanmean(koem_po4_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_po4_2016=nanstd(koem_po4_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_chl_2016=nanmean(koem_chl_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_chl_2016=nanstd(koem_chl_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_din_2016=nanmean(koem_din_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_din_2016=nanstd(koem_din_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_ss_2016=nanmean(koem_ss_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_ss_2016=nanstd(koem_ss_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_secchi_2016=nanmean(koem_secchi_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_secchi_2016=nanstd(koem_secchi_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_si_2016=nanmean(koem_si_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_si_2016=nanstd(koem_si_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_do_2016=nanmean(koem_do_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_do_2016=nanstd(koem_do_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% %bot
% d_b_t_2016=nanmean(koem_temp_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_t_2016=nanstd(koem_temp_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_s_2016=nanmean(koem_salt_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_s_2016=nanstd(koem_salt_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_nh4_2016=nanmean(koem_nh4_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_nh4_2016=nanstd(koem_nh4_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_no3_2016=nanmean(koem_no3_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_no3_2016=nanstd(koem_no3_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_po4_2016=nanmean(koem_po4_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_po4_2016=nanstd(koem_po4_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_chl_2016=nanmean(koem_chl_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_chl_2016=nanstd(koem_chl_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_din_2016=nanmean(koem_din_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_din_2016=nanstd(koem_din_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_ss_2016=nanmean(koem_ss_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_ss_2016=nanstd(koem_ss_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_si_2016=nanmean(koem_si_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_si_2016=nanstd(koem_si_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% d_b_do_2016=nanmean(koem_do_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% std_b_do_2016=nanstd(koem_do_b_raw(sp_9p_st,(j-1)*12+1:12*j),1);
% 

fig = figure; hold on;
        plot(find(isnan(d_t_2015)==0),d_t_2015(~isnan(d_t_2015)),'b','linew',2);
        plot(find(isnan(d_t_2016)==0),d_t_2016(~isnan(d_t_2016)),'r','linew',2);
        errorbar(d_t_2016,std_t_2016,'r','linew',2);
        errorbar(d_t_2015,std_t_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 temp']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp(^oC)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2016_KOEM_OBS_temp'),'-dpng') 
         
fig = figure; hold on;
        plot(find(isnan(d_b_t_2015)==0),d_b_t_2015(~isnan(d_t_2015)),'b','linew',2);
        plot(find(isnan(d_b_t_2016)==0),d_b_t_2016(~isnan(d_t_2016)),'r','linew',2);
        errorbar(d_b_t_2016,std_b_t_2016,'r','linew',2);
        errorbar(d_b_t_2015,std_b_t_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 temp bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp(^oC)','fontsize',13)
        grid on
                set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);
         set(gca,'fontsize',15)
         print(fig,strcat('2015vs2016_KOEM_OBS_temp_bot'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_s_2015)==0),d_s_2015(~isnan(d_s_2015)),'b','linew',2);
        plot(find(isnan(d_s_2016)==0),d_s_2016(~isnan(d_s_2016)),'r','linew',2);
        errorbar(d_s_2016,std_s_2016,'r','linew',2);
        errorbar(d_s_2015,std_s_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 salt']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt(psu)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
%         ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);      
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2016_KOEM_OBS_salt'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_b_s_2015)==0),d_b_s_2015(~isnan(d_b_s_2015)),'b','linew',2);
        plot(find(isnan(d_b_s_2016)==0),d_b_s_2016(~isnan(d_b_s_2016)),'r','linew',2);
        errorbar(d_b_s_2016,std_b_s_2016,'r','linew',2);
        errorbar(d_b_s_2015,std_b_s_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 salt bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt(psu)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 12])
%         ylim([0 1.2])
        set(gca,'xticklabel',1:12,'fontsize',10);        
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_salt_bot'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_nh4_2015)==0),d_nh4_2015(~isnan(d_nh4_2015)),'b','linew',2);
        plot(find(isnan(d_nh4_2016)==0),d_nh4_2016(~isnan(d_nh4_2016)),'r','linew',2);
        errorbar(d_nh4_2016,std_nh4_2016,'r','linew',2);
        errorbar(d_nh4_2015,std_nh4_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 8])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);    
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2016_KOEM_OBS_nh4'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_b_nh4_2015)==0),d_b_nh4_2015(~isnan(d_b_nh4_2015)),'b','linew',2);
        plot(find(isnan(d_b_nh4_2016)==0),d_b_nh4_2016(~isnan(d_b_nh4_2016)),'r','linew',2);
        errorbar(d_b_nh4_2016,std_b_nh4_2016,'r','linew',2);
        errorbar(d_b_nh4_2015,std_b_nh4_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 nh4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 8])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);        
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_nh4_bot'),'-dpng') 
        
        
fig = figure; hold on;
        plot(find(isnan(d_no3_2015)==0),d_no3_2015(~isnan(d_no3_2015)),'b','linew',2);
        plot(find(isnan(d_no3_2016)==0),d_no3_2016(~isnan(d_no3_2016)),'r','linew',2);
        errorbar(d_no3_2016,std_no3_2016,'r','linew',2);
        errorbar(d_no3_2015,std_no3_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 14])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_no3'),'-dpng');
        
fig = figure; hold on;
        plot(find(isnan(d_b_no3_2015)==0),d_b_no3_2015(~isnan(d_b_no3_2015)),'b','linew',2);
        plot(find(isnan(d_b_no3_2016)==0),d_b_no3_2016(~isnan(d_b_no3_2016)),'r','linew',2);
        errorbar(d_b_no3_2016,std_b_no3_2016,'r','linew',2);
        errorbar(d_b_no3_2015,std_b_no3_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 no3 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 14])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_no3_bot'),'-dpng')
        
        
fig = figure; hold on;
        plot(find(isnan(d_din_2015)==0),d_din_2015(~isnan(d_din_2015)),'b','linew',2);
        plot(find(isnan(d_din_2016)==0),d_din_2016(~isnan(d_din_2016)),'r','linew',2);
        errorbar(d_din_2016,std_din_2016,'r','linew',2);
        errorbar(d_din_2015,std_din_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 DIN']);
        xlabel('time(month)','fontsize',13)
        ylabel('DIN(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 23])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2016_KOEM_OBS_din'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_din_2015)==0),d_b_din_2015(~isnan(d_b_din_2015)),'b','linew',2);
        plot(find(isnan(d_b_din_2016)==0),d_b_din_2016(~isnan(d_b_din_2016)),'r','linew',2);
        errorbar(d_b_din_2016,std_b_din_2016,'r','linew',2);
        errorbar(d_b_din_2015,std_b_din_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 DIN bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('DIN(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 23])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);   
        print(fig,strcat('2015vs2016_KOEM_OBS_din_bot'),'-dpng')
        
 fig = figure; hold on;
        plot(find(isnan(d_si_2015)==0),d_si_2015(~isnan(d_si_2015)),'b','linew',2);
        plot(find(isnan(d_si_2016)==0),d_si_2016(~isnan(d_si_2016)),'r','linew',2);
        errorbar(d_si_2016,std_si_2016,'r','linew',2);
        errorbar(d_si_2015,std_si_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 si']);
        xlabel('time(month)','fontsize',13)
        ylabel('si(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 33])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_si'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_si_2015)==0),d_b_si_2015(~isnan(d_b_si_2015)),'b','linew',2);
        plot(find(isnan(d_b_si_2016)==0),d_b_si_2016(~isnan(d_b_si_2016)),'r','linew',2);
        errorbar(d_b_si_2016,std_b_si_2016,'r','linew',2);
        errorbar(d_b_si_2015,std_b_si_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 si bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('si(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 33])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_si_bot'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_po4_2015)==0),d_po4_2015(~isnan(d_po4_2015)),'b','linew',2);
        plot(find(isnan(d_po4_2016)==0),d_po4_2016(~isnan(d_po4_2016)),'r','linew',2);
        errorbar(d_po4_2016,std_po4_2016,'r','linew',2);
        errorbar(d_po4_2015,std_po4_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 po4']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 1.4])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_po4'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_po4_2015)==0),d_b_po4_2015(~isnan(d_b_po4_2015)),'b','linew',2);
        plot(find(isnan(d_b_po4_2016)==0),d_b_po4_2016(~isnan(d_b_po4_2016)),'r','linew',2);
        errorbar(d_b_po4_2016,std_b_po4_2016,'r','linew',2);
        errorbar(d_b_po4_2015,std_b_po4_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 po4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 1.4])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_po4_bot'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_do_2015)==0),d_do_2015(~isnan(d_do_2015)),'b','linew',2);
        plot(find(isnan(d_do_2016)==0),d_do_2016(~isnan(d_do_2016)),'r','linew',2);
        errorbar(d_do_2016,std_do_2016,'r','linew',2);
        errorbar(d_do_2015,std_do_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 DO']);
        xlabel('time(month)','fontsize',13)
        ylabel('DO (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_do'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_do_2015)==0),d_b_do_2015(~isnan(d_b_do_2015)),'b','linew',2);
        plot(find(isnan(d_b_do_2016)==0),d_b_do_2016(~isnan(d_b_do_2016)),'r','linew',2);
        errorbar(d_b_do_2016,std_b_do_2016,'r','linew',2);
        errorbar(d_b_do_2015,std_b_do_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 DO bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('DO (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
         print(fig,strcat('2015vs2016_KOEM_OBS_do_bot'),'-dpng')
        
        
 fig = figure; hold on;
        plot(find(isnan(d_chl_2015)==0),d_chl_2015(~isnan(d_chl_2015)),'b','linew',2);
        plot(find(isnan(d_chl_2016)==0),d_chl_2016(~isnan(d_chl_2016)),'r','linew',2);
        errorbar(d_chl_2016,std_chl_2016,'r','linew',2);
        errorbar(d_chl_2015,std_chl_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_chl'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_chl_2015)==0),d_b_chl_2015(~isnan(d_b_chl_2015)),'b','linew',2);
        plot(find(isnan(d_b_chl_2016)==0),d_b_chl_2016(~isnan(d_b_chl_2016)),'r','linew',2);
        errorbar(d_b_chl_2016,std_b_chl_2016,'r','linew',2);
        errorbar(d_b_chl_2015,std_b_chl_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 chl bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_chl_bot'),'-dpng')
        

fig = figure; hold on;
        plot(find(isnan(d_secchi_2015)==0),d_secchi_2015(~isnan(d_secchi_2015)),'b','linew',2);
        plot(find(isnan(d_secchi_2016)==0),d_secchi_2016(~isnan(d_secchi_2016)),'r','linew',2);
        errorbar(d_secchi_2016,std_secchi_2016,'r','linew',2);
        errorbar(d_secchi_2015,std_secchi_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 secchi']);
        xlabel('time(month)','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 5])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2016_KOEM_OBS_secchi'),'-dpng')
        
        
fig = figure; hold on;
        plot(find(isnan(d_ss_2015)==0),d_ss_2015(~isnan(d_ss_2015)),'b','linew',2);
        plot(find(isnan(d_ss_2016)==0),d_ss_2016(~isnan(d_ss_2016)),'r','linew',2);
        errorbar(d_ss_2016,std_ss_2016,'r','linew',2);
        errorbar(d_ss_2015,std_ss_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 ss']);
        xlabel('time(month)','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 39])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_ss'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_ss_2015)==0),d_b_ss_2015(~isnan(d_b_ss_2015)),'b','linew',2);
        plot(find(isnan(d_b_ss_2016)==0),d_b_ss_2016(~isnan(d_b_ss_2016)),'r','linew',2);
        errorbar(d_b_ss_2016,std_b_ss_2016,'r','linew',2);
        errorbar(d_b_ss_2015,std_b_ss_2015,'b','linew',2);
        title(['GY 2015 vs. 2016 ss bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 39])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_KOEM_OBS_ss_bot'),'-dpng')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4years (2008, 2013, 2015 vs. 2016)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fig = figure; hold on;
        plot(find(isnan(d_t_2015)==0),d_t_2015(~isnan(d_t_2015)),'b','linew',2);
        plot(find(isnan(d_t_2016)==0),d_t_2016(~isnan(d_t_2016)),'r','linew',2);
        errorbar(d_t_2016,std_t_2016,'r','linew',2);
        errorbar(d_t_2015,std_t_2015,'b','linew',2);
        plot(find(isnan(d_t_2008)==0),d_t_2008(~isnan(d_t_2008)),'k','linew',2);
        plot(find(isnan(d_t_2013)==0),d_t_2013(~isnan(d_t_2013)),'c','linew',2);
        errorbar(d_t_2008,std_t_2008,'k','linew',2);
        errorbar(d_t_2013,std_t_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 temp']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp(^oC)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);
        set(gca,'fontsize',15)
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_temp'),'-dpng') 
         
fig = figure; hold on;
        plot(find(isnan(d_b_t_2015)==0),d_b_t_2015(~isnan(d_t_2015)),'b','linew',2);
        plot(find(isnan(d_b_t_2016)==0),d_b_t_2016(~isnan(d_t_2016)),'r','linew',2);
        errorbar(d_b_t_2016,std_b_t_2016,'r','linew',2);
        errorbar(d_b_t_2015,std_b_t_2015,'b','linew',2);
        plot(find(isnan(d_b_t_2008)==0),d_b_t_2008(~isnan(d_b_t_2008)),'k','linew',2);
        plot(find(isnan(d_b_t_2013)==0),d_b_t_2013(~isnan(d_b_t_2013)),'c','linew',2);
        errorbar(d_b_t_2008,std_b_t_2008,'k','linew',2);
        errorbar(d_b_t_2013,std_b_t_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 temp bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp(^oC)','fontsize',13)
        grid on
                set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);
         set(gca,'fontsize',15)
         print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_temp_bot'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_s_2015)==0),d_s_2015(~isnan(d_s_2015)),'b','linew',2);
        plot(find(isnan(d_s_2016)==0),d_s_2016(~isnan(d_s_2016)),'r','linew',2);
        errorbar(d_s_2016,std_s_2016,'r','linew',2);
        errorbar(d_s_2015,std_s_2015,'b','linew',2);
        plot(find(isnan(d_s_2008)==0),d_s_2008(~isnan(d_s_2008)),'k','linew',2);
        plot(find(isnan(d_s_2013)==0),d_s_2013(~isnan(d_s_2013)),'c','linew',2);
        errorbar(d_s_2008,std_s_2008,'k','linew',2);
        errorbar(d_s_2013,std_s_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 salt']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt(psu)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
%         ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);      
        set(gca,'fontsize',15)
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_salt'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_b_s_2015)==0),d_b_s_2015(~isnan(d_b_s_2015)),'b','linew',2);
        plot(find(isnan(d_b_s_2016)==0),d_b_s_2016(~isnan(d_b_s_2016)),'r','linew',2);
        errorbar(d_b_s_2016,std_b_s_2016,'r','linew',2);
        errorbar(d_b_s_2015,std_b_s_2015,'b','linew',2);
        plot(find(isnan(d_b_s_2008)==0),d_b_s_2008(~isnan(d_b_s_2008)),'k','linew',2);
        plot(find(isnan(d_b_s_2013)==0),d_b_s_2013(~isnan(d_b_s_2013)),'c','linew',2);
        errorbar(d_b_s_2008,std_b_s_2008,'k','linew',2);
        errorbar(d_b_s_2013,std_b_s_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 salt bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt(psu)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 12])
%         ylim([0 1.2])
        set(gca,'xticklabel',1:12,'fontsize',10);        
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_salt_bot'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_nh4_2015)==0),d_nh4_2015(~isnan(d_nh4_2015)),'b','linew',2);
        plot(find(isnan(d_nh4_2016)==0),d_nh4_2016(~isnan(d_nh4_2016)),'r','linew',2);
        errorbar(d_nh4_2016,std_nh4_2016,'r','linew',2);
        errorbar(d_nh4_2015,std_nh4_2015,'b','linew',2);
        plot(find(isnan(d_nh4_2008)==0),d_nh4_2008(~isnan(d_nh4_2008)),'k','linew',2);
        plot(find(isnan(d_nh4_2013)==0),d_nh4_2013(~isnan(d_nh4_2013)),'c','linew',2);
        errorbar(d_nh4_2008,std_nh4_2008,'k','linew',2);
        errorbar(d_nh4_2013,std_nh4_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 8])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);    
        set(gca,'fontsize',15)
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_nh4'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_b_nh4_2015)==0),d_b_nh4_2015(~isnan(d_b_nh4_2015)),'b','linew',2);
        plot(find(isnan(d_b_nh4_2016)==0),d_b_nh4_2016(~isnan(d_b_nh4_2016)),'r','linew',2);
        errorbar(d_b_nh4_2016,std_b_nh4_2016,'r','linew',2);
        errorbar(d_b_nh4_2015,std_b_nh4_2015,'b','linew',2);
        plot(find(isnan(d_b_nh4_2008)==0),d_b_nh4_2008(~isnan(d_b_nh4_2008)),'k','linew',2);
        plot(find(isnan(d_b_nh4_2013)==0),d_b_nh4_2013(~isnan(d_b_nh4_2013)),'c','linew',2);
        errorbar(d_b_nh4_2008,std_b_nh4_2008,'k','linew',2);
        errorbar(d_b_nh4_2013,std_b_nh4_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 nh4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 8])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);        
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_nh4_bot'),'-dpng') 
        
        
fig = figure; hold on;
        plot(find(isnan(d_no3_2015)==0),d_no3_2015(~isnan(d_no3_2015)),'b','linew',2);
        plot(find(isnan(d_no3_2016)==0),d_no3_2016(~isnan(d_no3_2016)),'r','linew',2);
        errorbar(d_no3_2016,std_no3_2016,'r','linew',2);
        errorbar(d_no3_2015,std_no3_2015,'b','linew',2);
        plot(find(isnan(d_no3_2008)==0),d_no3_2008(~isnan(d_no3_2008)),'k','linew',2);
        plot(find(isnan(d_no3_2013)==0),d_no3_2013(~isnan(d_no3_2013)),'c','linew',2);
        errorbar(d_no3_2008,std_no3_2008,'k','linew',2);
        errorbar(d_no3_2013,std_no3_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 14])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_no3'),'-dpng');
        
fig = figure; hold on;
        plot(find(isnan(d_b_no3_2015)==0),d_b_no3_2015(~isnan(d_b_no3_2015)),'b','linew',2);
        plot(find(isnan(d_b_no3_2016)==0),d_b_no3_2016(~isnan(d_b_no3_2016)),'r','linew',2);
        errorbar(d_b_no3_2016,std_b_no3_2016,'r','linew',2);
        errorbar(d_b_no3_2015,std_b_no3_2015,'b','linew',2);
        plot(find(isnan(d_b_no3_2008)==0),d_b_no3_2008(~isnan(d_b_no3_2008)),'k','linew',2);
        plot(find(isnan(d_b_no3_2013)==0),d_b_no3_2013(~isnan(d_b_no3_2013)),'c','linew',2);
        errorbar(d_b_no3_2008,std_b_no3_2008,'k','linew',2);
        errorbar(d_b_no3_2013,std_b_no3_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 no3 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 14])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_no3_bot'),'-dpng')
        
        
fig = figure; hold on;
        plot(find(isnan(d_din_2015)==0),d_din_2015(~isnan(d_din_2015)),'b','linew',2);
        plot(find(isnan(d_din_2016)==0),d_din_2016(~isnan(d_din_2016)),'r','linew',2);
        errorbar(d_din_2016,std_din_2016,'r','linew',2);
        errorbar(d_din_2015,std_din_2015,'b','linew',2);
        plot(find(isnan(d_din_2008)==0),d_din_2008(~isnan(d_din_2008)),'k','linew',2);
        plot(find(isnan(d_din_2013)==0),d_din_2013(~isnan(d_din_2013)),'c','linew',2);
        errorbar(d_din_2008,std_din_2008,'k','linew',2);
        errorbar(d_din_2013,std_din_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 DIN']);
        xlabel('time(month)','fontsize',13)
        ylabel('DIN(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 23])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15)
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_din'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_din_2015)==0),d_b_din_2015(~isnan(d_b_din_2015)),'b','linew',2);
        plot(find(isnan(d_b_din_2016)==0),d_b_din_2016(~isnan(d_b_din_2016)),'r','linew',2);
        errorbar(d_b_din_2016,std_b_din_2016,'r','linew',2);
        errorbar(d_b_din_2015,std_b_din_2015,'b','linew',2);
        plot(find(isnan(d_b_din_2008)==0),d_b_din_2008(~isnan(d_b_din_2008)),'k','linew',2);
        plot(find(isnan(d_b_din_2013)==0),d_b_din_2013(~isnan(d_b_din_2013)),'c','linew',2);
        errorbar(d_b_din_2008,std_b_din_2008,'k','linew',2);
        errorbar(d_b_din_2013,std_b_din_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 DIN bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('DIN(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 23])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);   
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_din_bot'),'-dpng')
        
 fig = figure; hold on;
        plot(find(isnan(d_si_2015)==0),d_si_2015(~isnan(d_si_2015)),'b','linew',2);
        plot(find(isnan(d_si_2016)==0),d_si_2016(~isnan(d_si_2016)),'r','linew',2);
        errorbar(d_si_2016,std_si_2016,'r','linew',2);
        errorbar(d_si_2015,std_si_2015,'b','linew',2);
        plot(find(isnan(d_si_2008)==0),d_si_2008(~isnan(d_si_2008)),'k','linew',2);
        plot(find(isnan(d_si_2013)==0),d_si_2013(~isnan(d_si_2013)),'c','linew',2);
        errorbar(d_si_2008,std_si_2008,'k','linew',2);
        errorbar(d_si_2013,std_si_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 si']);
        xlabel('time(month)','fontsize',13)
        ylabel('si(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 33])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_si'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_si_2015)==0),d_b_si_2015(~isnan(d_b_si_2015)),'b','linew',2);
        plot(find(isnan(d_b_si_2016)==0),d_b_si_2016(~isnan(d_b_si_2016)),'r','linew',2);
        errorbar(d_b_si_2016,std_b_si_2016,'r','linew',2);
        errorbar(d_b_si_2015,std_b_si_2015,'b','linew',2);
        plot(find(isnan(d_b_si_2008)==0),d_b_si_2008(~isnan(d_b_si_2008)),'k','linew',2);
        plot(find(isnan(d_b_si_2013)==0),d_b_si_2013(~isnan(d_b_si_2013)),'c','linew',2);
        errorbar(d_b_si_2008,std_b_si_2008,'k','linew',2);
        errorbar(d_b_si_2013,std_b_si_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 si bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('si(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 33])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_si_bot'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_po4_2015)==0),d_po4_2015(~isnan(d_po4_2015)),'b','linew',2);
        plot(find(isnan(d_po4_2016)==0),d_po4_2016(~isnan(d_po4_2016)),'r','linew',2);
        errorbar(d_po4_2016,std_po4_2016,'r','linew',2);
        errorbar(d_po4_2015,std_po4_2015,'b','linew',2);
        plot(find(isnan(d_po4_2008)==0),d_po4_2008(~isnan(d_po4_2008)),'k','linew',2);
        plot(find(isnan(d_po4_2013)==0),d_po4_2013(~isnan(d_po4_2013)),'c','linew',2);
        errorbar(d_po4_2008,std_po4_2008,'k','linew',2);
        errorbar(d_po4_2013,std_po4_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 po4']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 1.4])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_po4'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_po4_2015)==0),d_b_po4_2015(~isnan(d_b_po4_2015)),'b','linew',2);
        plot(find(isnan(d_b_po4_2016)==0),d_b_po4_2016(~isnan(d_b_po4_2016)),'r','linew',2);
        errorbar(d_b_po4_2016,std_b_po4_2016,'r','linew',2);
        errorbar(d_b_po4_2015,std_b_po4_2015,'b','linew',2);
        plot(find(isnan(d_b_po4_2008)==0),d_b_po4_2008(~isnan(d_b_po4_2008)),'k','linew',2);
        plot(find(isnan(d_b_po4_2013)==0),d_b_po4_2013(~isnan(d_b_po4_2013)),'c','linew',2);
        errorbar(d_b_po4_2008,std_b_po4_2008,'k','linew',2);
        errorbar(d_b_po4_2013,std_b_po4_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 po4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 1.4])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_po4_bot'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_do_2015)==0),d_do_2015(~isnan(d_do_2015)),'b','linew',2);
        plot(find(isnan(d_do_2016)==0),d_do_2016(~isnan(d_do_2016)),'r','linew',2);
        errorbar(d_do_2016,std_do_2016,'r','linew',2);
        errorbar(d_do_2015,std_do_2015,'b','linew',2);
        plot(find(isnan(d_do_2008)==0),d_do_2008(~isnan(d_do_2008)),'k','linew',2);
        plot(find(isnan(d_do_2013)==0),d_do_2013(~isnan(d_do_2013)),'c','linew',2);
        errorbar(d_do_2008,std_do_2008,'k','linew',2);
        errorbar(d_do_2013,std_do_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 DO']);
        xlabel('time(month)','fontsize',13)
        ylabel('DO (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_do'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_do_2015)==0),d_b_do_2015(~isnan(d_b_do_2015)),'b','linew',2);
        plot(find(isnan(d_b_do_2016)==0),d_b_do_2016(~isnan(d_b_do_2016)),'r','linew',2);
        errorbar(d_b_do_2016,std_b_do_2016,'r','linew',2);
        errorbar(d_b_do_2015,std_b_do_2015,'b','linew',2);
        plot(find(isnan(d_b_do_2008)==0),d_b_do_2008(~isnan(d_b_do_2008)),'k','linew',2);
        plot(find(isnan(d_b_do_2013)==0),d_b_do_2013(~isnan(d_b_do_2013)),'c','linew',2);
        errorbar(d_b_do_2008,std_b_do_2008,'k','linew',2);
        errorbar(d_b_do_2013,std_b_do_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 DO bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('DO (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
         print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_do_bot'),'-dpng')
        
        
 fig = figure; hold on;
        plot(find(isnan(d_chl_2015)==0),d_chl_2015(~isnan(d_chl_2015)),'b','linew',2);
        plot(find(isnan(d_chl_2016)==0),d_chl_2016(~isnan(d_chl_2016)),'r','linew',2);
        errorbar(d_chl_2016,std_chl_2016,'r','linew',2);
        errorbar(d_chl_2015,std_chl_2015,'b','linew',2);
        plot(find(isnan(d_chl_2008)==0),d_chl_2008(~isnan(d_chl_2008)),'k','linew',2);
        plot(find(isnan(d_chl_2013)==0),d_chl_2013(~isnan(d_chl_2013)),'c','linew',2);
        errorbar(d_chl_2008,std_chl_2008,'k','linew',2);
        errorbar(d_chl_2013,std_chl_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 21])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_chl'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_chl_2015)==0),d_b_chl_2015(~isnan(d_b_chl_2015)),'b','linew',2);
        plot(find(isnan(d_b_chl_2016)==0),d_b_chl_2016(~isnan(d_b_chl_2016)),'r','linew',2);
        errorbar(d_b_chl_2016,std_b_chl_2016,'r','linew',2);
        errorbar(d_b_chl_2015,std_b_chl_2015,'b','linew',2);
        plot(find(isnan(d_b_chl_2008)==0),d_b_chl_2008(~isnan(d_b_chl_2008)),'k','linew',2);
        plot(find(isnan(d_b_chl_2013)==0),d_b_chl_2013(~isnan(d_b_chl_2013)),'c','linew',2);
        errorbar(d_b_chl_2008,std_b_chl_2008,'k','linew',2);
        errorbar(d_b_chl_2013,std_b_chl_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 chl bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 21])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_chl_bot'),'-dpng')
        

fig = figure; hold on;
        plot(find(isnan(d_secchi_2015)==0),d_secchi_2015(~isnan(d_secchi_2015)),'b','linew',2);
        plot(find(isnan(d_secchi_2016)==0),d_secchi_2016(~isnan(d_secchi_2016)),'r','linew',2);
        errorbar(d_secchi_2016,std_secchi_2016,'r','linew',2);
        errorbar(d_secchi_2015,std_secchi_2015,'b','linew',2);
        plot(find(isnan(d_secchi_2008)==0),d_secchi_2008(~isnan(d_secchi_2008)),'k','linew',2);
        plot(find(isnan(d_secchi_2013)==0),d_secchi_2013(~isnan(d_secchi_2013)),'c','linew',2);
        errorbar(d_secchi_2008,std_secchi_2008,'k','linew',2);
        errorbar(d_secchi_2013,std_secchi_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 secchi']);
        xlabel('time(month)','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 6])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15)
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_secchi'),'-dpng')
        
        
fig = figure; hold on;
        plot(find(isnan(d_ss_2015)==0),d_ss_2015(~isnan(d_ss_2015)),'b','linew',2);
        plot(find(isnan(d_ss_2016)==0),d_ss_2016(~isnan(d_ss_2016)),'r','linew',2);
        errorbar(d_ss_2016,std_ss_2016,'r','linew',2);
        errorbar(d_ss_2015,std_ss_2015,'b','linew',2);
        plot(find(isnan(d_ss_2008)==0),d_ss_2008(~isnan(d_ss_2008)),'k','linew',2);
        plot(find(isnan(d_ss_2013)==0),d_ss_2013(~isnan(d_ss_2013)),'c','linew',2);
        errorbar(d_ss_2008,std_ss_2008,'k','linew',2);
        errorbar(d_ss_2013,std_ss_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 ss']);
        xlabel('time(month)','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 39])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_ss'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_ss_2015)==0),d_b_ss_2015(~isnan(d_b_ss_2015)),'b','linew',2);
        plot(find(isnan(d_b_ss_2016)==0),d_b_ss_2016(~isnan(d_b_ss_2016)),'r','linew',2);
        errorbar(d_b_ss_2016,std_b_ss_2016,'r','linew',2);
        errorbar(d_b_ss_2015,std_b_ss_2015,'b','linew',2);
        plot(find(isnan(d_b_ss_2008)==0),d_b_ss_2008(~isnan(d_b_ss_2008)),'k','linew',2);
        plot(find(isnan(d_b_ss_2013)==0),d_b_ss_2013(~isnan(d_b_ss_2013)),'c','linew',2);
        errorbar(d_b_ss_2008,std_b_ss_2008,'k','linew',2);
        errorbar(d_b_ss_2013,std_b_ss_2013,'c','linew',2);
        title(['GY 2008, 2013, 2015 vs. 2016 ss bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 39])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015_2016_KOEM_OBS_ss_bot'),'-dpng')
        
        
        
        
        
%% load nut. 
clearvars plt_nh4* plt_no3* plt_po4*
i=27
plt_dis = (monthly_trans((i-1)*12+1 -2:12*i) .* 24*3600);
plt_nh4=sj.monthly_nh4((i-1)*12+1 -2:12*i) ./1000 .*14.006720.* (monthly_trans((i-1)*12+1 -2:12*i) .* 24*3600) ./ 10^6; 
plt_no3=sj.monthly_no3((i-1)*12+1 -2:12*i) ./1000 .*14.006720.* (monthly_trans((i-1)*12+1 -2:12*i) * 24*3600) ./ 10^6;
plt_po4=sj.monthly_po4((i-1)*12+1 -2:12*i) ./1000 .*30.973762.* (monthly_trans((i-1)*12+1 -2:12*i) * 24*3600) ./ 10^6;
j=28
plt_dis_16 = (monthly_trans((j-1)*12+1 -2:12*j) .* 24*3600);
plt_nh4_16=sj.monthly_nh4((j-1)*12+1 -2:12*j) ./1000 .*14.006720.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;
plt_no3_16=sj.monthly_no3((j-1)*12+1 -2:12*j) ./1000 .*14.006720.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;
plt_po4_16=sj.monthly_po4((j-1)*12+1 -2:12*j) ./1000 .*30.973762.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;

j=20
plt_dis_08 = (monthly_trans((j-1)*12+1 -2:12*j) .* 24*3600);
plt_nh4_08=sj.monthly_nh4((j-1)*12+1 -2:12*j) ./1000 .*14.006720.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;
plt_no3_08=sj.monthly_no3((j-1)*12+1 -2:12*j) ./1000 .*14.006720.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;
plt_po4_08=sj.monthly_po4((j-1)*12+1 -2:12*j) ./1000 .*30.973762.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;

j=25
plt_dis_13 = (monthly_trans((j-1)*12+1 -2:12*j) .* 24*3600);
plt_nh4_13=sj.monthly_nh4((j-1)*12+1 -2:12*j) ./1000 .*14.006720.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;
plt_no3_13=sj.monthly_no3((j-1)*12+1 -2:12*j) ./1000 .*14.006720.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;
plt_po4_13=sj.monthly_po4((j-1)*12+1 -2:12*j) ./1000 .*30.973762.* (monthly_trans((j-1)*12+1 -2:12*j).* 24*3600) ./ 10^6;


% plt_po4=sj.monthly_po4 ./1000 .*30.973762.* monthly_trans .* 31536000 ./ 10^6
fig = figure; hold on;
        plot(plt_no3,'bo','linew',2);
        plot(plt_no3,'b','linew',2);
        plot(plt_no3_16,'ro','linew',2);
        plot(plt_no3_16,'r','linew',2);
        title(['NO3-N load songjung']);
        xlabel('time(month)','fontsize',13)
        ylabel('load (ton/day)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
%         set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 16])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_river_no3_load'),'-dpng')
 
fig = figure; hold on;
        plot(plt_nh4,'bo','linew',2);
        plot(plt_nh4,'b','linew',2);
        plot(plt_nh4_16,'ro','linew',2);
        plot(plt_nh4_16,'r','linew',2);
        title(['NH4-N load songjung']);
        xlabel('time(month)','fontsize',13)
        ylabel('load (ton/day)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
%         set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 .7])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_river_nh4_load'),'-dpng')
        
               
fig = figure; hold on;
        plot(plt_po4,'bo','linew',2);
        plot(plt_po4,'b','linew',2);
        plot(plt_po4_16,'ro','linew',2);
        plot(plt_po4_16,'r','linew',2);
        title(['PO4-P load songjung']);
        xlabel('time(month)','fontsize',13)
        ylabel('load (ton/day)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
%         set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 .7])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2016_river_po4_load'),'-dpng')
 
 %% 3years river       
 
 fig = figure; hold on;
        plot(plt_dis ./ (24*3600),'bo','linew',2);
        plot(plt_dis ./ (24*3600),'b','linew',2);
        plot(plt_dis_16 ./ (24*3600),'ro','linew',2);
        plot(plt_dis_16 ./ (24*3600),'r','linew',2);
        plot(plt_dis_08 ./ (24*3600),'ko','linew',2);
        plot(plt_dis_08 ./ (24*3600),'k','linew',2);
        plot(plt_dis_13 ./ (24*3600),'co','linew',2);
        plot(plt_dis_13 ./ (24*3600),'c','linew',2);
        title(['river discharge songjung']);
        xlabel('time(month)','fontsize',13)
        ylabel('discharge (m^3/s)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:10:170);
        xlim([1 14])
        ylim([0 170])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015vs2016_river_discharge'),'-dpng')
 
        
 fig = figure; hold on;
        plot(plt_no3,'bo','linew',2);
        plot(plt_no3,'b','linew',2);
        plot(plt_no3_16,'ro','linew',2);
        plot(plt_no3_16,'r','linew',2);
        plot(plt_no3_08,'ko','linew',2);
        plot(plt_no3_08,'k','linew',2);
        plot(plt_no3_13,'co','linew',2);
        plot(plt_no3_13,'c','linew',2);
        title(['NO3-N load songjung']);
        xlabel('time(month)','fontsize',13)
        ylabel('load (ton/day)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:2:16);
        xlim([1 14])
        ylim([0 16])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015vs2016_river_no3_load'),'-dpng')
 
fig = figure; hold on;
        plot(plt_nh4,'bo','linew',2);
        plot(plt_nh4,'b','linew',2);
        plot(plt_nh4_16,'ro','linew',2);
        plot(plt_nh4_16,'r','linew',2);
        plot(plt_nh4_08,'ko','linew',2);
        plot(plt_nh4_08,'k','linew',2);
        plot(plt_nh4_13,'co','linew',2);
        plot(plt_nh4_13,'c','linew',2);
        title(['NH4-N load songjung']);
        xlabel('time(month)','fontsize',13)
        ylabel('load (ton/day)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
%         set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 .7])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015vs2016_river_nh4_load'),'-dpng')
        
               
fig = figure; hold on;
        plot(plt_po4,'bo','linew',2);
        plot(plt_po4,'b','linew',2);
        plot(plt_po4_16,'ro','linew',2);
        plot(plt_po4_16,'r','linew',2);
        plot(plt_po4_08,'ko','linew',2);
        plot(plt_po4_08,'k','linew',2);
        plot(plt_po4_13,'co','linew',2);
        plot(plt_po4_13,'c','linew',2);
        title(['PO4-P load songjung']);
        xlabel('time(month)','fontsize',13)
        ylabel('load (ton/day)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
%         set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 .7])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2008_2013_2015vs2016_river_po4_load'),'-dpng')       
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% 2015 vs. 2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fig = figure; hold on;
        plot(find(isnan(d_t_2015)==0),d_t_2015(~isnan(d_t_2015)),'b','linew',2);
        plot(find(isnan(d_t_2017)==0),d_t_2017(~isnan(d_t_2017)),'r','linew',2);
        errorbar(d_t_2017,std_t_2017,'r','linew',2);
        errorbar(d_t_2015,std_t_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 temp']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp(^oC)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2017_KOEM_OBS_temp'),'-dpng') 
         
fig = figure; hold on;
        plot(find(isnan(d_b_t_2015)==0),d_b_t_2015(~isnan(d_t_2015)),'b','linew',2);
        plot(find(isnan(d_b_t_2017)==0),d_b_t_2017(~isnan(d_t_2017)),'r','linew',2);
        errorbar(d_b_t_2017,std_b_t_2017,'r','linew',2);
        errorbar(d_b_t_2015,std_b_t_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 temp bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp(^oC)','fontsize',13)
        grid on
                set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);
         set(gca,'fontsize',15)
         print(fig,strcat('2015vs2017_KOEM_OBS_temp_bot'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_s_2015)==0),d_s_2015(~isnan(d_s_2015)),'b','linew',2);
        plot(find(isnan(d_s_2017)==0),d_s_2017(~isnan(d_s_2017)),'r','linew',2);
        errorbar(d_s_2017,std_s_2017,'r','linew',2);
        errorbar(d_s_2015,std_s_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 salt']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt(psu)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
%         ylim([5 30])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);      
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2017_KOEM_OBS_salt'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_b_s_2015)==0),d_b_s_2015(~isnan(d_b_s_2015)),'b','linew',2);
        plot(find(isnan(d_b_s_2017)==0),d_b_s_2017(~isnan(d_b_s_2017)),'r','linew',2);
        errorbar(d_b_s_2017,std_b_s_2017,'r','linew',2);
        errorbar(d_b_s_2015,std_b_s_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 salt bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt(psu)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 12])
%         ylim([0 1.2])
        set(gca,'xticklabel',1:12,'fontsize',10);        
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_salt_bot'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_nh4_2015)==0),d_nh4_2015(~isnan(d_nh4_2015)),'b','linew',2);
        plot(find(isnan(d_nh4_2017)==0),d_nh4_2017(~isnan(d_nh4_2017)),'r','linew',2);
        errorbar(d_nh4_2017,std_nh4_2017,'r','linew',2);
        errorbar(d_nh4_2015,std_nh4_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 8])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);    
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2017_KOEM_OBS_nh4'),'-dpng') 
        
fig = figure; hold on;
        plot(find(isnan(d_b_nh4_2015)==0),d_b_nh4_2015(~isnan(d_b_nh4_2015)),'b','linew',2);
        plot(find(isnan(d_b_nh4_2017)==0),d_b_nh4_2017(~isnan(d_b_nh4_2017)),'r','linew',2);
        errorbar(d_b_nh4_2017,std_b_nh4_2017,'r','linew',2);
        errorbar(d_b_nh4_2015,std_b_nh4_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 nh4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 8])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);        
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_nh4_bot'),'-dpng') 
        
        
fig = figure; hold on;
        plot(find(isnan(d_no3_2015)==0),d_no3_2015(~isnan(d_no3_2015)),'b','linew',2);
        plot(find(isnan(d_no3_2017)==0),d_no3_2017(~isnan(d_no3_2017)),'r','linew',2);
        errorbar(d_no3_2017,std_no3_2017,'r','linew',2);
        errorbar(d_no3_2015,std_no3_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 no3']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 14])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_no3'),'-dpng');
        
fig = figure; hold on;
        plot(find(isnan(d_b_no3_2015)==0),d_b_no3_2015(~isnan(d_b_no3_2015)),'b','linew',2);
        plot(find(isnan(d_b_no3_2017)==0),d_b_no3_2017(~isnan(d_b_no3_2017)),'r','linew',2);
        errorbar(d_b_no3_2017,std_b_no3_2017,'r','linew',2);
        errorbar(d_b_no3_2015,std_b_no3_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 no3 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 14])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_no3_bot'),'-dpng')
        
        
fig = figure; hold on;
        plot(find(isnan(d_din_2015)==0),d_din_2015(~isnan(d_din_2015)),'b','linew',2);
        plot(find(isnan(d_din_2017)==0),d_din_2017(~isnan(d_din_2017)),'r','linew',2);
        errorbar(d_din_2017,std_din_2017,'r','linew',2);
        errorbar(d_din_2015,std_din_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 DIN']);
        xlabel('time(month)','fontsize',13)
        ylabel('DIN(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 23])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2017_KOEM_OBS_din'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_din_2015)==0),d_b_din_2015(~isnan(d_b_din_2015)),'b','linew',2);
        plot(find(isnan(d_b_din_2017)==0),d_b_din_2017(~isnan(d_b_din_2017)),'r','linew',2);
        errorbar(d_b_din_2017,std_b_din_2017,'r','linew',2);
        errorbar(d_b_din_2015,std_b_din_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 DIN bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('DIN(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 23])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);   
        print(fig,strcat('2015vs2017_KOEM_OBS_din_bot'),'-dpng')
        
 fig = figure; hold on;
        plot(find(isnan(d_si_2015)==0),d_si_2015(~isnan(d_si_2015)),'b','linew',2);
        plot(find(isnan(d_si_2017)==0),d_si_2017(~isnan(d_si_2017)),'r','linew',2);
        errorbar(d_si_2017,std_si_2017,'r','linew',2);
        errorbar(d_si_2015,std_si_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 si']);
        xlabel('time(month)','fontsize',13)
        ylabel('si(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 33])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_si'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_si_2015)==0),d_b_si_2015(~isnan(d_b_si_2015)),'b','linew',2);
        plot(find(isnan(d_b_si_2017)==0),d_b_si_2017(~isnan(d_b_si_2017)),'r','linew',2);
        errorbar(d_b_si_2017,std_b_si_2017,'r','linew',2);
        errorbar(d_b_si_2015,std_b_si_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 si bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('si(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 33])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_si_bot'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_po4_2015)==0),d_po4_2015(~isnan(d_po4_2015)),'b','linew',2);
        plot(find(isnan(d_po4_2017)==0),d_po4_2017(~isnan(d_po4_2017)),'r','linew',2);
        errorbar(d_po4_2017,std_po4_2017,'r','linew',2);
        errorbar(d_po4_2015,std_po4_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 po4']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 1.4])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_po4'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_po4_2015)==0),d_b_po4_2015(~isnan(d_b_po4_2015)),'b','linew',2);
        plot(find(isnan(d_b_po4_2017)==0),d_b_po4_2017(~isnan(d_b_po4_2017)),'r','linew',2);
        errorbar(d_b_po4_2017,std_b_po4_2017,'r','linew',2);
        errorbar(d_b_po4_2015,std_b_po4_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 po4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4(umol/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 1.4])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_po4_bot'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_do_2015)==0),d_do_2015(~isnan(d_do_2015)),'b','linew',2);
        plot(find(isnan(d_do_2017)==0),d_do_2017(~isnan(d_do_2017)),'r','linew',2);
        errorbar(d_do_2017,std_do_2017,'r','linew',2);
        errorbar(d_do_2015,std_do_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 DO']);
        xlabel('time(month)','fontsize',13)
        ylabel('DO (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_do'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_do_2015)==0),d_b_do_2015(~isnan(d_b_do_2015)),'b','linew',2);
        plot(find(isnan(d_b_do_2017)==0),d_b_do_2017(~isnan(d_b_do_2017)),'r','linew',2);
        errorbar(d_b_do_2017,std_b_do_2017,'r','linew',2);
        errorbar(d_b_do_2015,std_b_do_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 DO bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('DO (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
         print(fig,strcat('2015vs2017_KOEM_OBS_do_bot'),'-dpng')
        
        
 fig = figure; hold on;
        plot(find(isnan(d_chl_2015)==0),d_chl_2015(~isnan(d_chl_2015)),'b','linew',2);
        plot(find(isnan(d_chl_2017)==0),d_chl_2017(~isnan(d_chl_2017)),'r','linew',2);
        errorbar(d_chl_2017,std_chl_2017,'r','linew',2);
        errorbar(d_chl_2015,std_chl_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_chl'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_chl_2015)==0),d_b_chl_2015(~isnan(d_b_chl_2015)),'b','linew',2);
        plot(find(isnan(d_b_chl_2017)==0),d_b_chl_2017(~isnan(d_b_chl_2017)),'r','linew',2);
        errorbar(d_b_chl_2017,std_b_chl_2017,'r','linew',2);
        errorbar(d_b_chl_2015,std_b_chl_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 chl bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 12])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_chl_bot'),'-dpng')
        

fig = figure; hold on;
        plot(find(isnan(d_secchi_2015)==0),d_secchi_2015(~isnan(d_secchi_2015)),'b','linew',2);
        plot(find(isnan(d_secchi_2017)==0),d_secchi_2017(~isnan(d_secchi_2017)),'r','linew',2);
        errorbar(d_secchi_2017,std_secchi_2017,'r','linew',2);
        errorbar(d_secchi_2015,std_secchi_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 secchi']);
        xlabel('time(month)','fontsize',13)
        ylabel('secchi (m)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        xlim([1 14])
        ylim([0 5])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15)
        print(fig,strcat('2015vs2017_KOEM_OBS_secchi'),'-dpng')
        
        
fig = figure; hold on;
        plot(find(isnan(d_ss_2015)==0),d_ss_2015(~isnan(d_ss_2015)),'b','linew',2);
        plot(find(isnan(d_ss_2017)==0),d_ss_2017(~isnan(d_ss_2017)),'r','linew',2);
        errorbar(d_ss_2017,std_ss_2017,'r','linew',2);
        errorbar(d_ss_2015,std_ss_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 ss']);
        xlabel('time(month)','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 39])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);       
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_ss'),'-dpng')
        
fig = figure; hold on;
        plot(find(isnan(d_b_ss_2015)==0),d_b_ss_2015(~isnan(d_b_ss_2015)),'b','linew',2);
        plot(find(isnan(d_b_ss_2017)==0),d_b_ss_2017(~isnan(d_b_ss_2017)),'r','linew',2);
        errorbar(d_b_ss_2017,std_b_ss_2017,'r','linew',2);
        errorbar(d_b_ss_2015,std_b_ss_2015,'b','linew',2);
        title(['GY 2015 vs. 2017 ss bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('ss (mg/L)','fontsize',13)
        grid on
        set(gca,'xtick',1:31);
        set(gca,'ytick',0:5:40);
        xlim([1 14])
        ylim([0 39])
        set(gca,'xticklabel',[11:12,1:12],'fontsize',10);         
        set(gca,'fontsize',15);
        print(fig,strcat('2015vs2017_KOEM_OBS_ss_bot'),'-dpng')
        
