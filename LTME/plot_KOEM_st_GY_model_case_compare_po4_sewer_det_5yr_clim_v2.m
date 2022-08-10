close all; clear;clc

case_num = [1:3];
sigsig=[2;3;];
% 1: spatial & time std.
% 2: time std.
% 3: spatial std.
 dir_name ={'all_std','t_std','sp_std'};


cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
ltme =load('LTME_observation_data_v2.mat'); %LTME obs.

% cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_2001_koem_daily\spmean
% load('koem_result_processed_std_spmean.mat'); % KOEM obs.
% % mpdata= load('mpdata_result.mat');

cd D:\장기생태\Dynamic\result\1997
c1997 = load('1997_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\1998
c1998 = load('1998_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\1999
c1999 = load('1999_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2000
c2000 = load('2000_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2001
c2001 = load('2001_mp_p_sewer_det_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2002
c2002 = load('2002_mp_p_sewer_det_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2003
c2003 = load('2003_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2004
c2004 = load('2004_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2005
c2005 = load('2005_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2006
c2006 = load('2006_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2007
c2007 = load('2007_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2008
c2008 = load('2008_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

cd D:\장기생태\Dynamic\result\2009
c2009 = load('2009_mp_p_sewer_det_f_result.mat','obm_*','gy_*','jj_*','sgy_*','egy_*','obs_std_*','eom_d_each','sp_gy_*');

% first_5;
c1997_eom = [1 (c1997.eom_d_each(1:end-1) + 1)];
c1998_eom = [(c1997.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(1:end-1) + 1)];
c1999_eom = [(c1997.eom_d_each(end)+c1998.eom_d_each(end)+1) (c1997.eom_d_each(end)+c1998.eom_d_each(end)+c1999.eom_d_each(1:end-1) + 1)];
% second_5;
c2000_eom = [1 (c2000.eom_d_each(1:end-1)+ 1)];
c2001_eom = [(c2000.eom_d_each(end)+1) (c2000.eom_d_each(end)+c2001.eom_d_each(1:end-1)+ 1)];
c2002_eom = [(c2000.eom_d_each(end)+c2001.eom_d_each(end)+1) (c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(1:end-1) + 1)];
c2003_eom = [(c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(end)+1) (c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(1:end-1) + 1)];
c2004_eom = [(c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+1) (c2000.eom_d_each(end)+c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+c2004.eom_d_each(1:end-1) + 1)];
% third_5;
c2005_eom = [1 (c2005.eom_d_each(1:end-1)+ 1)];
c2006_eom = [(c2005.eom_d_each(end)+1) (c2005.eom_d_each(end)+c2006.eom_d_each(1:end-1)+ 1)];
c2007_eom = [(c2005.eom_d_each(end)+c2006.eom_d_each(end)+1) (c2005.eom_d_each(end)+c2006.eom_d_each(end)+c2007.eom_d_each(1:end-1) + 1)];
c2008_eom = [(c2005.eom_d_each(end)+c2006.eom_d_each(end)+c2007.eom_d_each(end)+1) (c2005.eom_d_each(end)+c2006.eom_d_each(end)+c2007.eom_d_each(end)+c2008.eom_d_each(1:end-1) + 1)];
c2009_eom = [(c2005.eom_d_each(end)+c2006.eom_d_each(end)+c2007.eom_d_each(end)+c2008.eom_d_each(end)+1) (c2005.eom_d_each(end)+c2006.eom_d_each(end)+c2007.eom_d_each(end)+c2008.eom_d_each(end)+c2009.eom_d_each(1:end-1) + 1)];

% c2005_eom = [(c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end)+c2004.eom_d_each(end)+1) (c2001.eom_d_each(end)+c2002.eom_d_each(end)+c2003.eom_d_each(end) + c2004.eom_d_each(end) + c2005.eom_d_each(1:end-1)+ 1)];


com_eom1 = [c1997_eom c1998_eom c1999_eom ];
com_eom2 = [c2000_eom c2001_eom c2002_eom c2003_eom c2004_eom];
com_eom3 = [c2005_eom c2006_eom c2007_eom c2008_eom c2009_eom];

% % for i = 1:length(t_tIc)+length(t_tIc2)
% %     if i < 13
% %         com_t_tic{i}=t_tIc{i}
% %     else
% %         com_t_tic{i}=t_tIc2{i-12}
% %     end
% % end

clearvars com_t_tic
% com_t_tic = {'2001.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2002.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2003.01','02','03','04','05','06','07','08','09','10','11','12',...}
%             '2004.01','02','03','04','05','06','07','08','09','10','11','12'}
com_t_tic = {'1997','1998','1999','2000','2001', ...
            '2002', ...
            '2003'};
        
c1997.sp_gy_no3(:,60)=[]; c1997.sp_gy_chl(:,60)=[]; c1997.sp_gy_nh4(:,60)=[]; c1997.sp_gy_temp(:,60)=[]; c1997.sp_gy_salt(:,60)=[];
c1998.sp_gy_no3(:,60)=[]; c1998.sp_gy_chl(:,60)=[]; c1998.sp_gy_nh4(:,60)=[]; c1998.sp_gy_temp(:,60)=[]; c1998.sp_gy_salt(:,60)=[];
c1999.sp_gy_no3(:,60)=[]; c1999.sp_gy_chl(:,60)=[]; c1999.sp_gy_nh4(:,60)=[]; c1999.sp_gy_temp(:,60)=[]; c1999.sp_gy_salt(:,60)=[];
c2000.sp_gy_no3(:,60)=[]; c2000.sp_gy_chl(:,60)=[]; c2000.sp_gy_nh4(:,60)=[]; c2000.sp_gy_temp(:,60)=[]; c2000.sp_gy_salt(:,60)=[];

c1997.sp_gy_no3_b(:,60)=[]; c1997.sp_gy_chl_b(:,60)=[]; c1997.sp_gy_nh4_b(:,60)=[]; c1997.sp_gy_temp_b(:,60)=[]; c1997.sp_gy_salt_b(:,60)=[];
c1998.sp_gy_no3_b(:,60)=[]; c1998.sp_gy_chl_b(:,60)=[]; c1998.sp_gy_nh4_b(:,60)=[]; c1998.sp_gy_temp_b(:,60)=[]; c1998.sp_gy_salt_b(:,60)=[];
c1999.sp_gy_no3_b(:,60)=[]; c1999.sp_gy_chl_b(:,60)=[]; c1999.sp_gy_nh4_b(:,60)=[]; c1999.sp_gy_temp_b(:,60)=[]; c1999.sp_gy_salt_b(:,60)=[];
c2000.sp_gy_no3_b(:,60)=[]; c2000.sp_gy_chl_b(:,60)=[]; c2000.sp_gy_nh4_b(:,60)=[]; c2000.sp_gy_temp_b(:,60)=[]; c2000.sp_gy_salt_b(:,60)=[];

c2004.sp_gy_no3(:,60)=[]; c2004.sp_gy_chl(:,60)=[]; c2004.sp_gy_nh4(:,60)=[]; c2004.sp_gy_temp(:,60)=[]; c2004.sp_gy_salt(:,60)=[];
c2005.sp_gy_no3(:,60)=[]; c2005.sp_gy_chl(:,60)=[]; c2005.sp_gy_nh4(:,60)=[]; c2005.sp_gy_temp(:,60)=[]; c2005.sp_gy_salt(:,60)=[];
c2006.sp_gy_no3(:,60)=[]; c2006.sp_gy_chl(:,60)=[]; c2006.sp_gy_nh4(:,60)=[]; c2006.sp_gy_temp(:,60)=[]; c2006.sp_gy_salt(:,60)=[];
c2007.sp_gy_no3(:,60)=[]; c2007.sp_gy_chl(:,60)=[]; c2007.sp_gy_nh4(:,60)=[]; c2007.sp_gy_temp(:,60)=[]; c2007.sp_gy_salt(:,60)=[];
c2008.sp_gy_no3(:,60)=[]; c2008.sp_gy_chl(:,60)=[]; c2008.sp_gy_nh4(:,60)=[]; c2008.sp_gy_temp(:,60)=[]; c2008.sp_gy_salt(:,60)=[];
c2009.sp_gy_no3(:,60)=[]; c2009.sp_gy_chl(:,60)=[]; c2009.sp_gy_nh4(:,60)=[]; c2009.sp_gy_temp(:,60)=[]; c2009.sp_gy_salt(:,60)=[];

c2004.sp_gy_no3_b(:,60)=[]; c2004.sp_gy_chl_b(:,60)=[]; c2004.sp_gy_nh4_b(:,60)=[]; c2004.sp_gy_temp_b(:,60)=[]; c2004.sp_gy_salt_b(:,60)=[];
c2005.sp_gy_no3_b(:,60)=[]; c2005.sp_gy_chl_b(:,60)=[]; c2005.sp_gy_nh4_b(:,60)=[]; c2005.sp_gy_temp_b(:,60)=[]; c2005.sp_gy_salt_b(:,60)=[];
c2006.sp_gy_no3_b(:,60)=[]; c2006.sp_gy_chl_b(:,60)=[]; c2006.sp_gy_nh4_b(:,60)=[]; c2006.sp_gy_temp_b(:,60)=[]; c2006.sp_gy_salt_b(:,60)=[];
c2007.sp_gy_no3_b(:,60)=[]; c2007.sp_gy_chl_b(:,60)=[]; c2007.sp_gy_nh4_b(:,60)=[]; c2007.sp_gy_temp_b(:,60)=[]; c2007.sp_gy_salt_b(:,60)=[];
c2008.sp_gy_no3_b(:,60)=[]; c2008.sp_gy_chl_b(:,60)=[]; c2008.sp_gy_nh4_b(:,60)=[]; c2008.sp_gy_temp_b(:,60)=[]; c2008.sp_gy_salt_b(:,60)=[];
c2009.sp_gy_no3_b(:,60)=[]; c2009.sp_gy_chl_b(:,60)=[]; c2009.sp_gy_nh4_b(:,60)=[]; c2009.sp_gy_temp_b(:,60)=[]; c2009.sp_gy_salt_b(:,60)=[];

c1997.gy_no3(60)=[]; c1997.gy_chl(60)=[]; c1997.gy_nh4(60)=[]; c1997.gy_temp(60)=[]; c1997.gy_salt(60)=[];
c1998.gy_no3(60)=[]; c1998.gy_chl(60)=[]; c1998.gy_nh4(60)=[]; c1998.gy_temp(60)=[]; c1998.gy_salt(60)=[];
c1999.gy_no3(60)=[]; c1999.gy_chl(60)=[]; c1999.gy_nh4(60)=[]; c1999.gy_temp(60)=[]; c1999.gy_salt(60)=[];
c2000.gy_no3(60)=[]; c2000.gy_chl(60)=[]; c2000.gy_nh4(60)=[]; c2000.gy_temp(60)=[]; c2000.gy_salt(60)=[];

c1997.gy_no3_b(60)=[]; c1997.gy_chl_b(60)=[]; c1997.gy_nh4_b(60)=[]; c1997.gy_temp_b(60)=[]; c1997.gy_salt_b(60)=[];
c1998.gy_no3_b(60)=[]; c1998.gy_chl_b(60)=[]; c1998.gy_nh4_b(60)=[]; c1998.gy_temp_b(60)=[]; c1998.gy_salt_b(60)=[];
c1999.gy_no3_b(60)=[]; c1999.gy_chl_b(60)=[]; c1999.gy_nh4_b(60)=[]; c1999.gy_temp_b(60)=[]; c1999.gy_salt_b(60)=[];
c2000.gy_no3_b(60)=[]; c2000.gy_chl_b(60)=[]; c2000.gy_nh4_b(60)=[]; c2000.gy_temp_b(60)=[]; c2000.gy_salt_b(60)=[];

c2004.gy_no3(60)=[]; c2004.gy_chl(60)=[]; c2004.gy_nh4(60)=[]; c2004.gy_temp(60)=[]; c2004.gy_salt(60)=[];
c2005.gy_no3(60)=[]; c2005.gy_chl(60)=[]; c2005.gy_nh4(60)=[]; c2005.gy_temp(60)=[]; c2005.gy_salt(60)=[];
c2006.gy_no3(60)=[]; c2006.gy_chl(60)=[]; c2006.gy_nh4(60)=[]; c2006.gy_temp(60)=[]; c2006.gy_salt(60)=[];
c2007.gy_no3(60)=[]; c2007.gy_chl(60)=[]; c2007.gy_nh4(60)=[]; c2007.gy_temp(60)=[]; c2007.gy_salt(60)=[];
c2008.gy_no3(60)=[]; c2008.gy_chl(60)=[]; c2008.gy_nh4(60)=[]; c2008.gy_temp(60)=[]; c2008.gy_salt(60)=[];
c2009.gy_no3(60)=[]; c2009.gy_chl(60)=[]; c2009.gy_nh4(60)=[]; c2009.gy_temp(60)=[]; c2009.gy_salt(60)=[];

c2004.gy_no3_b(60)=[]; c2004.gy_chl_b(60)=[]; c2004.gy_nh4_b(60)=[]; c2004.gy_temp_b(60)=[]; c2004.gy_salt_b(60)=[];
c2005.gy_no3_b(60)=[]; c2005.gy_chl_b(60)=[]; c2005.gy_nh4_b(60)=[]; c2005.gy_temp_b(60)=[]; c2005.gy_salt_b(60)=[];
c2006.gy_no3_b(60)=[]; c2006.gy_chl_b(60)=[]; c2006.gy_nh4_b(60)=[]; c2006.gy_temp_b(60)=[]; c2006.gy_salt_b(60)=[];
c2007.gy_no3_b(60)=[]; c2007.gy_chl_b(60)=[]; c2007.gy_nh4_b(60)=[]; c2007.gy_temp_b(60)=[]; c2007.gy_salt_b(60)=[];
c2008.gy_no3_b(60)=[]; c2008.gy_chl_b(60)=[]; c2008.gy_nh4_b(60)=[]; c2008.gy_temp_b(60)=[]; c2008.gy_salt_b(60)=[];
c2009.gy_no3_b(60)=[]; c2009.gy_chl_b(60)=[]; c2009.gy_nh4_b(60)=[]; c2009.gy_temp_b(60)=[]; c2009.gy_salt_b(60)=[];
        
%% 1st - 5yr        
mer_no3_1 = [c1997.gy_no3; c1998.gy_no3; c1999.gy_no3;];
mer_nh4_1 = [c1997.gy_nh4; c1998.gy_nh4; c1999.gy_nh4;];
mer_chl_1 = [c1997.gy_chl; c1998.gy_chl; c1999.gy_chl;];
mer_temp_1 = [c1997.gy_temp; c1998.gy_temp; c1999.gy_temp;];
mer_salt_1 = [c1997.gy_salt; c1998.gy_salt; c1999.gy_salt;];

%bot merged
mer_no3_b_1 = [c1997.gy_no3_b; c1998.gy_no3_b; c1999.gy_no3_b;];
mer_nh4_b_1 = [c1997.gy_nh4_b; c1998.gy_nh4_b; c1999.gy_nh4_b;];
mer_chl_b_1 = [c1997.gy_chl_b; c1998.gy_chl_b; c1999.gy_chl_b;];
mer_temp_b_1 = [c1997.gy_temp_b; c1998.gy_temp_b; c1999.gy_temp_b;];
mer_salt_b_1 = [c1997.gy_salt_b; c1998.gy_salt_b; c1999.gy_salt_b;];
        
% sptatial value
%surf merged
mer_sp_no3_1 = [c1997.sp_gy_no3 c1998.sp_gy_no3 c1999.sp_gy_no3 ];
mer_sp_nh4_1 = [c1997.sp_gy_nh4 c1998.sp_gy_nh4 c1999.sp_gy_nh4 ];
mer_sp_chl_1 = [c1997.sp_gy_chl c1998.sp_gy_chl c1999.sp_gy_chl];
mer_sp_temp_1 = [c1997.sp_gy_temp c1998.sp_gy_temp c1999.sp_gy_temp ];
mer_sp_salt_1 = [c1997.sp_gy_salt c1998.sp_gy_salt c1999.sp_gy_salt];
%bot merged
mer_sp_no3_b_1 = [c1997.sp_gy_no3_b c1998.sp_gy_no3_b c1999.sp_gy_no3_b ];
mer_sp_nh4_b_1 = [c1997.sp_gy_nh4_b c1998.sp_gy_nh4_b c1999.sp_gy_nh4_b ];
mer_sp_chl_b_1 = [c1997.sp_gy_chl_b c1998.sp_gy_chl_b c1999.sp_gy_chl_b ];
mer_sp_temp_b_1 = [c1997.sp_gy_temp_b c1998.sp_gy_temp_b c1999.sp_gy_temp_b ];
mer_sp_salt_b_1 = [c1997.sp_gy_salt_b c1998.sp_gy_salt_b c1999.sp_gy_salt_b ];

%% 2nd - 5yr   
mer_no3_2 = [c2000.gy_no3; c2001.gy_no3; c2002.gy_no3; c2003.gy_no3; c2004.gy_no3;];
mer_nh4_2 = [c2000.gy_nh4; c2001.gy_nh4; c2002.gy_nh4; c2003.gy_nh4; c2004.gy_nh4;];
mer_chl_2 = [c2000.gy_chl; c2001.gy_chl; c2002.gy_chl; c2003.gy_chl; c2004.gy_chl;];
mer_temp_2 = [c2000.gy_temp; c2001.gy_temp; c2002.gy_temp; c2003.gy_temp; c2004.gy_temp;];
mer_salt_2 = [c2000.gy_salt; c2001.gy_salt; c2002.gy_salt; c2003.gy_salt; c2004.gy_salt;];

%bot merged
mer_no3_b_2 = [c2000.gy_no3_b; c2001.gy_no3_b; c2002.gy_no3_b; c2003.gy_no3_b; c2004.gy_no3_b;];
mer_nh4_b_2 = [c2000.gy_nh4_b; c2001.gy_nh4_b; c2002.gy_nh4_b; c2003.gy_nh4_b; c2004.gy_nh4_b;];
mer_chl_b_2 = [c2000.gy_chl_b; c2001.gy_chl_b; c2002.gy_chl_b; c2003.gy_chl_b; c2004.gy_chl_b;];
mer_temp_b_2 = [c2000.gy_temp_b; c2001.gy_temp_b; c2002.gy_temp_b; c2003.gy_temp_b; c2004.gy_temp_b;];
mer_salt_b_2 = [c2000.gy_salt_b; c2001.gy_salt_b; c2002.gy_salt_b; c2003.gy_salt_b; c2004.gy_salt_b;];

% sptatial value
%surf merged
mer_sp_no3_2 = [c2000.sp_gy_no3 c2001.sp_gy_no3 c2002.sp_gy_no3 c2003.sp_gy_no3 c2004.sp_gy_no3];
mer_sp_nh4_2 = [c2000.sp_gy_nh4 c2001.sp_gy_nh4 c2002.sp_gy_nh4 c2003.sp_gy_nh4 c2004.sp_gy_nh4];
mer_sp_chl_2 = [c2000.sp_gy_chl c2001.sp_gy_chl c2002.sp_gy_chl c2003.sp_gy_chl c2004.sp_gy_chl];
mer_sp_temp_2 = [c2000.sp_gy_temp c2001.sp_gy_temp c2002.sp_gy_temp c2003.sp_gy_temp c2004.sp_gy_temp];
mer_sp_salt_2 = [c2000.sp_gy_salt c2001.sp_gy_salt c2002.sp_gy_salt c2003.sp_gy_salt c2004.sp_gy_salt];
%bot merged
mer_sp_no3_b_2 = [c2000.sp_gy_no3_b c2001.sp_gy_no3_b c2002.sp_gy_no3_b c2003.sp_gy_no3_b c2004.sp_gy_no3_b];
mer_sp_nh4_b_2 = [c2000.sp_gy_nh4_b c2001.sp_gy_nh4_b c2002.sp_gy_nh4_b c2003.sp_gy_nh4_b c2004.sp_gy_nh4_b];
mer_sp_chl_b_2 = [c2000.sp_gy_chl_b c2001.sp_gy_chl_b c2002.sp_gy_chl_b c2003.sp_gy_chl_b c2004.sp_gy_chl_b];
mer_sp_temp_b_2 = [c2000.sp_gy_temp_b c2001.sp_gy_temp_b c2002.sp_gy_temp_b c2003.sp_gy_temp_b c2004.sp_gy_temp_b];
mer_sp_salt_b_2 = [c2000.sp_gy_salt_b c2001.sp_gy_salt_b c2002.sp_gy_salt_b c2003.sp_gy_salt_b c2004.sp_gy_salt_b];

%% 3rd - 5yr   
mer_no3_3 = [c2005.gy_no3; c2006.gy_no3; c2007.gy_no3; c2008.gy_no3; c2009.gy_no3;];
mer_nh4_3 = [c2005.gy_nh4; c2006.gy_nh4; c2007.gy_nh4; c2008.gy_nh4; c2009.gy_nh4;];
mer_chl_3 = [c2005.gy_chl; c2006.gy_chl; c2007.gy_chl; c2008.gy_chl; c2009.gy_chl;];
mer_temp_3 = [c2005.gy_temp; c2006.gy_temp; c2007.gy_temp; c2008.gy_temp; c2009.gy_temp;];
mer_salt_3 = [c2005.gy_salt; c2006.gy_salt; c2007.gy_salt; c2008.gy_salt; c2009.gy_salt;];

%bot merged
mer_no3_b_3 = [c2005.gy_no3_b; c2006.gy_no3_b; c2007.gy_no3_b; c2008.gy_no3_b; c2009.gy_no3_b;];
mer_nh4_b_3 = [c2005.gy_nh4_b; c2006.gy_nh4_b; c2007.gy_nh4_b; c2008.gy_nh4_b; c2009.gy_nh4_b;];
mer_chl_b_3 = [c2005.gy_chl_b; c2006.gy_chl_b; c2007.gy_chl_b; c2008.gy_chl_b; c2009.gy_chl_b;];
mer_temp_b_3 = [c2005.gy_temp_b; c2006.gy_temp_b; c2007.gy_temp_b; c2008.gy_temp_b; c2009.gy_temp_b;];
mer_salt_b_3 = [c2005.gy_salt_b; c2006.gy_salt_b; c2007.gy_salt_b; c2008.gy_salt_b; c2009.gy_salt_b;];

% sptatial value
%surf merged
mer_sp_no3_3 = [c2005.sp_gy_no3 c2006.sp_gy_no3 c2007.sp_gy_no3 c2008.sp_gy_no3 c2009.sp_gy_no3];
mer_sp_nh4_3 = [c2005.sp_gy_nh4 c2006.sp_gy_nh4 c2007.sp_gy_nh4 c2008.sp_gy_nh4 c2009.sp_gy_nh4];
mer_sp_chl_3 = [c2005.sp_gy_chl c2006.sp_gy_chl c2007.sp_gy_chl c2008.sp_gy_chl c2009.sp_gy_chl];
mer_sp_temp_3 = [c2005.sp_gy_temp c2006.sp_gy_temp c2007.sp_gy_temp c2008.sp_gy_temp c2009.sp_gy_temp];
mer_sp_salt_3 = [c2005.sp_gy_salt c2006.sp_gy_salt c2007.sp_gy_salt c2008.sp_gy_salt c2009.sp_gy_salt];
%bot merged
mer_sp_no3_b_3 = [c2005.sp_gy_no3_b c2006.sp_gy_no3_b c2007.sp_gy_no3_b c2008.sp_gy_no3_b c2009.sp_gy_no3_b];
mer_sp_nh4_b_3 = [c2005.sp_gy_nh4_b c2006.sp_gy_nh4_b c2007.sp_gy_nh4_b c2008.sp_gy_nh4_b c2009.sp_gy_nh4_b];
mer_sp_chl_b_3 = [c2005.sp_gy_chl_b c2006.sp_gy_chl_b c2007.sp_gy_chl_b c2008.sp_gy_chl_b c2009.sp_gy_chl_b];
mer_sp_temp_b_3 = [c2005.sp_gy_temp_b c2006.sp_gy_temp_b c2007.sp_gy_temp_b c2008.sp_gy_temp_b c2009.sp_gy_temp_b];
mer_sp_salt_b_3 = [c2005.sp_gy_salt_b c2006.sp_gy_salt_b c2007.sp_gy_salt_b c2008.sp_gy_salt_b c2009.sp_gy_salt_b];


for case_num = 1:3
    for ixx=1:2
        clearvars -except ixx case_num sigsig c199* c200* ltme *_eom dir_name mer_*
        
        sig=sigsig(ixx)

cd D:\장기생태\Dynamic\koem
% clim1=load('koem_climate_3regime_v2_2sig_gy_to_2003_12.mat'); % from 'plot_KOEM_st_GY_model_compare_horizontal_plt_climate_3regime_v2_gy_2sig.m'
% clim1=load('koem_climate_3regime_v2_to_2003_12.mat');  % from 'plot_KOEM_st_GY_model_compare_horizontal_plt_seasonal_climate_3regime_v2.m'
% clim1=load('koem_climate_3regime_v2_3sig_gy_to_2003_12.mat');  % from 'plot_KOEM_st_seasonal_climate_3regime_v3_spatial_sigma_extract.m'

! clim1=load(['koem_climate_3regime_v2_',num2str(sig),'sig_gy_to_2003_12.mat']);

% 1997~2018 : 22yr
% 1st regime : 1997 ~ 2003.12 (from koem nh4 shift)
% 2nd regime : 2004.01 ~ 2009.12 (from koem no3 shift)
% 3rd regime : 2010.01 ~ 2018.12 (from koem no3 shift)

% 5yr cutting
% 1997~1999 : first_5;
% 2000~2004 : second_5;
% 2005~2009 : third_5;

% % mpdata.gy_temp

%make time tick
for i = 1:12
t_tIc{i} = {num2str(i,'%02d')};
if i == 1
    t_tIc{i} = {['2001.' num2str(i,'%02d')]};
end
end

for i = 1:12
t_tIc2{i} = {num2str(i,'%02d')};
if i == 1
    t_tIc2{i} = {['2002.' num2str(i,'%02d')]};
end
end

% for i = 1:12
% t_tIc2{i} = {num2str(i,'%02d')};
% if i == 1
%     t_tIc3{i} = {['2003.' num2str(i,'%02d')]};
% end
% end



for i =  1:365
%% model 
% time & spatial std
    %% 1st
    clearvars size2
    size2=size(mer_sp_no3_1(:,i:365:end),2);
    regm_std_no3_1(i)=std(reshape(mer_sp_no3_1(:,i:365:end),1,9*size2));
    regm_std_nh4_1(i)=std(reshape(mer_sp_nh4_1(:,i:365:end),1,9*size2));
    regm_std_chl_1(i)=std(reshape(mer_sp_chl_1(:,i:365:end),1,9*size2));
    regm_std_temp_1(i)=std(reshape(mer_sp_temp_1(:,i:365:end),1,9*size2));
    regm_std_salt_1(i)=std(reshape(mer_sp_salt_1(:,i:365:end),1,9*size2));
  
    clearvars size2
    size2=size(mer_sp_no3_1(:,i:365:end),2);
    regm_std_no3_b_1(i)=std(reshape(mer_sp_no3_b_1(:,i:365:end),1,9*size2));
    regm_std_nh4_b_1(i)=std(reshape(mer_sp_nh4_b_1(:,i:365:end),1,9*size2));
    regm_std_chl_b_1(i)=std(reshape(mer_sp_chl_b_1(:,i:365:end),1,9*size2));
    regm_std_temp_b_1(i)=std(reshape(mer_sp_temp_b_1(:,i:365:end),1,9*size2));
    regm_std_salt_b_1(i)=std(reshape(mer_sp_salt_b_1(:,i:365:end),1,9*size2));

    %% 2nd
    clearvars size2
    size2=size(mer_sp_no3_2(:,i:365:end),2);
    regm_std_no3_2(i)=std(reshape(mer_sp_no3_2(:,i:365:end),1,9*size2));
    regm_std_nh4_2(i)=std(reshape(mer_sp_nh4_2(:,i:365:end),1,9*size2));
    regm_std_chl_2(i)=std(reshape(mer_sp_chl_2(:,i:365:end),1,9*size2));
    regm_std_temp_2(i)=std(reshape(mer_sp_temp_2(:,i:365:end),1,9*size2));
    regm_std_salt_2(i)=std(reshape(mer_sp_salt_2(:,i:365:end),1,9*size2));
  

    clearvars size2
    size2=size(mer_sp_no3_2(:,i:365:end),2);
    regm_std_no3_b_2(i)=std(reshape(mer_sp_no3_b_2(:,i:365:end),1,9*size2));
    regm_std_nh4_b_2(i)=std(reshape(mer_sp_nh4_b_2(:,i:365:end),1,9*size2));
    regm_std_chl_b_2(i)=std(reshape(mer_sp_chl_b_2(:,i:365:end),1,9*size2));
    regm_std_temp_b_2(i)=std(reshape(mer_sp_temp_b_2(:,i:365:end),1,9*size2));
    regm_std_salt_b_2(i)=std(reshape(mer_sp_salt_b_2(:,i:365:end),1,9*size2));
    
    %% 3rd
    clearvars size2
    size2=size(mer_sp_no3_3(:,i:365:end),2);
    regm_std_no3_3(i)=std(reshape(mer_sp_no3_3(:,i:365:end),1,9*size2));
    regm_std_nh4_3(i)=std(reshape(mer_sp_nh4_3(:,i:365:end),1,9*size2));
    regm_std_chl_3(i)=std(reshape(mer_sp_chl_3(:,i:365:end),1,9*size2));
    regm_std_temp_3(i)=std(reshape(mer_sp_temp_3(:,i:365:end),1,9*size2));
    regm_std_salt_3(i)=std(reshape(mer_sp_salt_3(:,i:365:end),1,9*size2));
  

    clearvars size2
    size2=size(mer_sp_no3_3(:,i:365:end),2);
    regm_std_no3_b_3(i)=std(reshape(mer_sp_no3_b_3(:,i:365:end),1,9*size2));
    regm_std_nh4_b_3(i)=std(reshape(mer_sp_nh4_b_3(:,i:365:end),1,9*size2));
    regm_std_chl_b_3(i)=std(reshape(mer_sp_chl_b_3(:,i:365:end),1,9*size2));
    regm_std_temp_b_3(i)=std(reshape(mer_sp_temp_b_3(:,i:365:end),1,9*size2));
    regm_std_salt_b_3(i)=std(reshape(mer_sp_salt_b_3(:,i:365:end),1,9*size2));
    
% time std  
    %% 1st
    spm_regm_std_no3_1(i)=std(nanmean(mer_sp_no3_1(:,i:365:end),1));
    spm_regm_std_nh4_1(i)=std(nanmean(mer_sp_nh4_1(:,i:365:end),1));
    spm_regm_std_chl_1(i)=std(nanmean(mer_sp_chl_1(:,i:365:end),1));
    spm_regm_std_temp_1(i)=std(nanmean(mer_sp_temp_1(:,i:365:end),1));
    spm_regm_std_salt_1(i)=std(nanmean(mer_sp_salt_1(:,i:365:end),1));
    
    spm_regm_std_no3_b_1(i)=std(nanmean(mer_sp_no3_b_1(:,i:365:end),1));
    spm_regm_std_nh4_b_1(i)=std(nanmean(mer_sp_nh4_b_1(:,i:365:end),1));
    spm_regm_std_chl_b_1(i)=std(nanmean(mer_sp_chl_b_1(:,i:365:end),1));
    spm_regm_std_temp_b_1(i)=std(nanmean(mer_sp_temp_b_1(:,i:365:end),1));
    spm_regm_std_salt_b_1(i)=std(nanmean(mer_sp_salt_b_1(:,i:365:end),1));
    %% 2nd
    spm_regm_std_no3_2(i)=std(nanmean(mer_sp_no3_2(:,i:365:end),1));
    spm_regm_std_nh4_2(i)=std(nanmean(mer_sp_nh4_2(:,i:365:end),1));
    spm_regm_std_chl_2(i)=std(nanmean(mer_sp_chl_2(:,i:365:end),1));
    spm_regm_std_temp_2(i)=std(nanmean(mer_sp_temp_2(:,i:365:end),1));
    spm_regm_std_salt_2(i)=std(nanmean(mer_sp_salt_2(:,i:365:end),1));
    
    spm_regm_std_no3_b_2(i)=std(nanmean(mer_sp_no3_b_2(:,i:365:end),1));
    spm_regm_std_nh4_b_2(i)=std(nanmean(mer_sp_nh4_b_2(:,i:365:end),1));
    spm_regm_std_chl_b_2(i)=std(nanmean(mer_sp_chl_b_2(:,i:365:end),1));
    spm_regm_std_temp_b_2(i)=std(nanmean(mer_sp_temp_b_2(:,i:365:end),1));
    spm_regm_std_salt_b_2(i)=std(nanmean(mer_sp_salt_b_2(:,i:365:end),1));
    %% 3rd
    spm_regm_std_no3_3(i)=std(nanmean(mer_sp_no3_3(:,i:365:end),1));
    spm_regm_std_nh4_3(i)=std(nanmean(mer_sp_nh4_3(:,i:365:end),1));
    spm_regm_std_chl_3(i)=std(nanmean(mer_sp_chl_3(:,i:365:end),1));
    spm_regm_std_temp_3(i)=std(nanmean(mer_sp_temp_3(:,i:365:end),1));
    spm_regm_std_salt_3(i)=std(nanmean(mer_sp_salt_3(:,i:365:end),1));
    
    spm_regm_std_no3_b_3(i)=std(nanmean(mer_sp_no3_b_3(:,i:365:end),1));
    spm_regm_std_nh4_b_3(i)=std(nanmean(mer_sp_nh4_b_3(:,i:365:end),1));
    spm_regm_std_chl_b_3(i)=std(nanmean(mer_sp_chl_b_3(:,i:365:end),1));
    spm_regm_std_temp_b_3(i)=std(nanmean(mer_sp_temp_b_3(:,i:365:end),1));
    spm_regm_std_salt_b_3(i)=std(nanmean(mer_sp_salt_b_3(:,i:365:end),1));
    
%spatial std  
    %% 1st
    tm_regm_std_no3_1(i)=std(nanmean(mer_sp_no3_1(:,i:365:end),2));
    tm_regm_std_nh4_1(i)=std(nanmean(mer_sp_nh4_1(:,i:365:end),2));
    tm_regm_std_chl_1(i)=std(nanmean(mer_sp_chl_1(:,i:365:end),2));
    tm_regm_std_temp_1(i)=std(nanmean(mer_sp_temp_1(:,i:365:end),2));
    tm_regm_std_salt_1(i)=std(nanmean(mer_sp_salt_1(:,i:365:end),2));
    
    tm_regm_std_no3_b_1(i)=std(nanmean(mer_sp_no3_b_1(:,i:365:end),2));
    tm_regm_std_nh4_b_1(i)=std(nanmean(mer_sp_nh4_b_1(:,i:365:end),2));
    tm_regm_std_chl_b_1(i)=std(nanmean(mer_sp_chl_b_1(:,i:365:end),2));
    tm_regm_std_temp_b_1(i)=std(nanmean(mer_sp_temp_b_1(:,i:365:end),2));
    tm_regm_std_salt_b_1(i)=std(nanmean(mer_sp_salt_b_1(:,i:365:end),2));
    %% 2nd
    tm_regm_std_no3_2(i)=std(nanmean(mer_sp_no3_2(:,i:365:end),2));
    tm_regm_std_nh4_2(i)=std(nanmean(mer_sp_nh4_2(:,i:365:end),2));
    tm_regm_std_chl_2(i)=std(nanmean(mer_sp_chl_2(:,i:365:end),2));
    tm_regm_std_temp_2(i)=std(nanmean(mer_sp_temp_2(:,i:365:end),2));
    tm_regm_std_salt_2(i)=std(nanmean(mer_sp_salt_2(:,i:365:end),2));
    
    tm_regm_std_no3_b_2(i)=std(nanmean(mer_sp_no3_b_2(:,i:365:end),2));
    tm_regm_std_nh4_b_2(i)=std(nanmean(mer_sp_nh4_b_2(:,i:365:end),2));
    tm_regm_std_chl_b_2(i)=std(nanmean(mer_sp_chl_b_2(:,i:365:end),2));
    tm_regm_std_temp_b_2(i)=std(nanmean(mer_sp_temp_b_2(:,i:365:end),2));
    tm_regm_std_salt_b_2(i)=std(nanmean(mer_sp_salt_b_2(:,i:365:end),2));
    %% 3rd
    tm_regm_std_no3_3(i)=std(nanmean(mer_sp_no3_3(:,i:365:end),2));
    tm_regm_std_nh4_3(i)=std(nanmean(mer_sp_nh4_3(:,i:365:end),2));
    tm_regm_std_chl_3(i)=std(nanmean(mer_sp_chl_3(:,i:365:end),2));
    tm_regm_std_temp_3(i)=std(nanmean(mer_sp_temp_3(:,i:365:end),2));
    tm_regm_std_salt_3(i)=std(nanmean(mer_sp_salt_3(:,i:365:end),2));
    
    tm_regm_std_no3_b_3(i)=std(nanmean(mer_sp_no3_b_3(:,i:365:end),2));
    tm_regm_std_nh4_b_3(i)=std(nanmean(mer_sp_nh4_b_3(:,i:365:end),2));
    tm_regm_std_chl_b_3(i)=std(nanmean(mer_sp_chl_b_3(:,i:365:end),2));
    tm_regm_std_temp_b_3(i)=std(nanmean(mer_sp_temp_b_3(:,i:365:end),2));
    tm_regm_std_salt_b_3(i)=std(nanmean(mer_sp_salt_b_3(:,i:365:end),2));  
    
end   

return 

% %  dir_name ={'all_std','t_std','sp_std'};
cd(['D:\장기생태\Dynamic\KOEM\compare_sewer\1st_5yr\1st_5yr_',num2str(sig),'sig_',dir_name{case_num}]);

if case_num == 1
fig = figure; hold on;
       