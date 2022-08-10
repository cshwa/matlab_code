close all; clear;clc

cd D:\장기생태\observation\관측\장기생태프로젝트\남해_DB
ltme =load('LTME_observation_data_v2.mat'); %LTME obs.

% cd D:\장기생태\Dynamic\KOEM\gy_2001\gy_1997~2010_koem_daily\spmean
% load('koem_result_processed_std_spmean.mat'); % KOEM obs.
% % mpdata= load('mpdata_result.mat');
cd D:\장기생태\Dynamic\result\1997
c1997 = load('1997_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\1998
c1998 = load('1998_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\1999
c1999 = load('1999_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2000
c2000 = load('2000_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2001
c2001 = load('2001_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2002
c2002 = load('2002_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2003
c2003 = load('2003_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2004
c2004 = load('2004_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2005
c2005 = load('2005_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2006
c2006 = load('2006_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2007
c2007 = load('2007_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2008
c2008 = load('2008_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2009
c2009 = load('2009_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');
cd D:\장기생태\Dynamic\result\2010
c2010 = load('2010_monthly_result.mat','gy_*','jj_*','sgy_*','egy_*','std_*','sp_gy_*');

cd D:\장기생태\Dynamic\result
% clim1=load('koem_compare_data_gy_2sig_monthly_spatial_mean.mat');
clim1=load('koem_compare_data_gy_3sig_monthly_spatial_mean.mat');

% cd D:\장기생태\Dynamic\KOEM\gy_2001\phymin_half_saturation
% phymin= load('phymin_result.mat');
% 
% cd D:\장기생태\Dynamic\KOEM\gy_2001\mp_p_ws_half_saturation
% mp_p_ws= load('mp_p_ws_result.mat');

% % mpdata.gy_temp
cd D:\장기생태\Dynamic\KOEM\compare_sewer\det_3clim_monthly_3sig


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
    t_tIc2{i} = {['1997~2010.' num2str(i,'%02d')]};
end
end

% for i = 1:12
% t_tIc2{i} = {num2str(i,'%02d')};
% if i == 1
%     t_tIc3{i} = {['2003.' num2str(i,'%02d')]};
% end
% end

%% variable
%surf merged
mer_sp_no3 = [c1997.sp_gy_no3 c1998.sp_gy_no3 c1999.sp_gy_no3 c2000.sp_gy_no3 c2001.sp_gy_no3 ...
    c2002.sp_gy_no3 c2003.sp_gy_no3 c2004.sp_gy_no3 c2005.sp_gy_no3 c2006.sp_gy_no3 c2007.sp_gy_no3 ...
    c2008.sp_gy_no3 c2009.sp_gy_no3 c2010.sp_gy_no3];
mer_sp_nh4 = [c1997.sp_gy_nh4 c1998.sp_gy_nh4 c1999.sp_gy_nh4 c2000.sp_gy_nh4 c2001.sp_gy_nh4 ...
    c2002.sp_gy_nh4 c2003.sp_gy_nh4 c2004.sp_gy_nh4 c2005.sp_gy_nh4 c2006.sp_gy_nh4 c2007.sp_gy_nh4 ...
    c2008.sp_gy_nh4 c2009.sp_gy_nh4 c2010.sp_gy_nh4];
mer_sp_chl = [c1997.sp_gy_chl c1998.sp_gy_chl c1999.sp_gy_chl c2000.sp_gy_chl c2001.sp_gy_chl ...
    c2002.sp_gy_chl c2003.sp_gy_chl c2004.sp_gy_chl c2005.sp_gy_chl c2006.sp_gy_chl c2007.sp_gy_chl ...
    c2008.sp_gy_chl c2009.sp_gy_chl c2010.sp_gy_chl];
mer_sp_temp = [c1997.sp_gy_temp c1998.sp_gy_temp c1999.sp_gy_temp c2000.sp_gy_temp c2001.sp_gy_temp ...
    c2002.sp_gy_temp c2003.sp_gy_temp c2004.sp_gy_temp c2005.sp_gy_temp c2006.sp_gy_temp c2007.sp_gy_temp ...
    c2008.sp_gy_temp c2009.sp_gy_temp c2010.sp_gy_temp];
mer_sp_salt = [c1997.sp_gy_salt c1998.sp_gy_salt c1999.sp_gy_salt c2000.sp_gy_salt c2001.sp_gy_salt ...
    c2002.sp_gy_salt c2003.sp_gy_salt c2004.sp_gy_salt c2005.sp_gy_salt c2006.sp_gy_salt c2007.sp_gy_salt ...
    c2008.sp_gy_salt c2009.sp_gy_salt c2010.sp_gy_salt];
%bot merged
mer_sp_no3_b = [c1997.sp_gy_no3_b c1998.sp_gy_no3_b c1999.sp_gy_no3_b c2000.sp_gy_no3_b c2001.sp_gy_no3_b ...
    c2002.sp_gy_no3_b c2003.sp_gy_no3_b c2004.sp_gy_no3_b c2005.sp_gy_no3_b c2006.sp_gy_no3_b c2007.sp_gy_no3_b ...
    c2008.sp_gy_no3_b c2009.sp_gy_no3_b c2010.sp_gy_no3_b];
mer_sp_nh4_b = [c1997.sp_gy_nh4_b c1998.sp_gy_nh4_b c1999.sp_gy_nh4_b c2000.sp_gy_nh4_b c2001.sp_gy_nh4_b ...
    c2002.sp_gy_nh4_b c2003.sp_gy_nh4_b c2004.sp_gy_nh4_b c2005.sp_gy_nh4_b c2006.sp_gy_nh4_b c2007.sp_gy_nh4_b ...
    c2008.sp_gy_nh4_b c2009.sp_gy_nh4_b c2010.sp_gy_nh4_b];
mer_sp_chl_b = [c1997.sp_gy_chl_b c1998.sp_gy_chl_b c1999.sp_gy_chl_b c2000.sp_gy_chl_b c2001.sp_gy_chl_b ...
    c2002.sp_gy_chl_b c2003.sp_gy_chl_b c2004.sp_gy_chl_b c2005.sp_gy_chl_b c2006.sp_gy_chl_b c2007.sp_gy_chl_b ...
    c2008.sp_gy_chl_b c2009.sp_gy_chl_b c2010.sp_gy_chl_b];
mer_sp_temp_b = [c1997.sp_gy_temp_b c1998.sp_gy_temp_b c1999.sp_gy_temp_b c2000.sp_gy_temp_b c2001.sp_gy_temp_b ...
    c2002.sp_gy_temp_b c2003.sp_gy_temp_b c2004.sp_gy_temp_b c2005.sp_gy_temp_b c2006.sp_gy_temp_b c2007.sp_gy_temp_b ...
    c2008.sp_gy_temp_b c2009.sp_gy_temp_b c2010.sp_gy_temp_b];
mer_sp_salt_b = [c1997.sp_gy_salt_b c1998.sp_gy_salt_b c1999.sp_gy_salt_b c2000.sp_gy_salt_b c2001.sp_gy_salt_b ...
    c2002.sp_gy_salt_b c2003.sp_gy_salt_b c2004.sp_gy_salt_b c2005.sp_gy_salt_b c2006.sp_gy_salt_b c2007.sp_gy_salt_b ...
    c2008.sp_gy_salt_b c2009.sp_gy_salt_b c2010.sp_gy_salt_b];

%surf merged
mer_no3 = [c1997.gy_no3 c1998.gy_no3 c1999.gy_no3 c2000.gy_no3 c2001.gy_no3 ...
    c2002.gy_no3 c2003.gy_no3 c2004.gy_no3 c2005.gy_no3 c2006.gy_no3 c2007.gy_no3 ...
    c2008.gy_no3 c2009.gy_no3 c2010.gy_no3];
mer_nh4 = [c1997.gy_nh4 c1998.gy_nh4 c1999.gy_nh4 c2000.gy_nh4 c2001.gy_nh4 ...
    c2002.gy_nh4 c2003.gy_nh4 c2004.gy_nh4 c2005.gy_nh4 c2006.gy_nh4 c2007.gy_nh4 ...
    c2008.gy_nh4 c2009.gy_nh4 c2010.gy_nh4];
mer_chl = [c1997.gy_chl c1998.gy_chl c1999.gy_chl c2000.gy_chl c2001.gy_chl ...
    c2002.gy_chl c2003.gy_chl c2004.gy_chl c2005.gy_chl c2006.gy_chl c2007.gy_chl ...
    c2008.gy_chl c2009.gy_chl c2010.gy_chl];
mer_temp = [c1997.gy_temp c1998.gy_temp c1999.gy_temp c2000.gy_temp c2001.gy_temp ...
    c2002.gy_temp c2003.gy_temp c2004.gy_temp c2005.gy_temp c2006.gy_temp c2007.gy_temp ...
    c2008.gy_temp c2009.gy_temp c2010.gy_temp];
mer_salt = [c1997.gy_salt c1998.gy_salt c1999.gy_salt c2000.gy_salt c2001.gy_salt ...
    c2002.gy_salt c2003.gy_salt c2004.gy_salt c2005.gy_salt c2006.gy_salt c2007.gy_salt ...
    c2008.gy_salt c2009.gy_salt c2010.gy_salt];
%bot merged
mer_no3_b = [c1997.gy_no3_b c1998.gy_no3_b c1999.gy_no3_b c2000.gy_no3_b c2001.gy_no3_b ...
    c2002.gy_no3_b c2003.gy_no3_b c2004.gy_no3_b c2005.gy_no3_b c2006.gy_no3_b c2007.gy_no3_b ...
    c2008.gy_no3_b c2009.gy_no3_b c2010.gy_no3_b];
mer_nh4_b = [c1997.gy_nh4_b c1998.gy_nh4_b c1999.gy_nh4_b c2000.gy_nh4_b c2001.gy_nh4_b ...
    c2002.gy_nh4_b c2003.gy_nh4_b c2004.gy_nh4_b c2005.gy_nh4_b c2006.gy_nh4_b c2007.gy_nh4_b ...
    c2008.gy_nh4_b c2009.gy_nh4_b c2010.gy_nh4_b];
mer_chl_b = [c1997.gy_chl_b c1998.gy_chl_b c1999.gy_chl_b c2000.gy_chl_b c2001.gy_chl_b ...
    c2002.gy_chl_b c2003.gy_chl_b c2004.gy_chl_b c2005.gy_chl_b c2006.gy_chl_b c2007.gy_chl_b ...
    c2008.gy_chl_b c2009.gy_chl_b c2010.gy_chl_b];
mer_temp_b = [c1997.gy_temp_b c1998.gy_temp_b c1999.gy_temp_b c2000.gy_temp_b c2001.gy_temp_b ...
    c2002.gy_temp_b c2003.gy_temp_b c2004.gy_temp_b c2005.gy_temp_b c2006.gy_temp_b c2007.gy_temp_b ...
    c2008.gy_temp_b c2009.gy_temp_b c2010.gy_temp_b];
mer_salt_b = [c1997.gy_salt_b c1998.gy_salt_b c1999.gy_salt_b c2000.gy_salt_b c2001.gy_salt_b ...
    c2002.gy_salt_b c2003.gy_salt_b c2004.gy_salt_b c2005.gy_salt_b c2006.gy_salt_b c2007.gy_salt_b ...
    c2008.gy_salt_b c2009.gy_salt_b c2010.gy_salt_b];

%% std _ model
%surf merged
mer_std_no3 = [c1997.std_gy_no3 c1998.std_gy_no3 c1999.std_gy_no3 c2000.std_gy_no3 c2001.std_gy_no3 ...
    c2002.std_gy_no3 c2003.std_gy_no3 c2004.std_gy_no3 c2005.std_gy_no3 c2006.std_gy_no3 c2007.std_gy_no3 ...
    c2008.std_gy_no3 c2009.std_gy_no3 c2010.std_gy_no3];
mer_std_nh4 = [c1997.std_gy_nh4 c1998.std_gy_nh4 c1999.std_gy_nh4 c2000.std_gy_nh4 c2001.std_gy_nh4 ...
    c2002.std_gy_nh4 c2003.std_gy_nh4 c2004.std_gy_nh4 c2005.std_gy_nh4 c2006.std_gy_nh4 c2007.std_gy_nh4 ...
    c2008.std_gy_nh4 c2009.std_gy_nh4 c2010.std_gy_nh4];
mer_std_chl = [c1997.std_gy_chl c1998.std_gy_chl c1999.std_gy_chl c2000.std_gy_chl c2001.std_gy_chl ...
    c2002.std_gy_chl c2003.std_gy_chl c2004.std_gy_chl c2005.std_gy_chl c2006.std_gy_chl c2007.std_gy_chl ...
    c2008.std_gy_chl c2009.std_gy_chl c2010.std_gy_chl];
mer_std_temp = [c1997.std_gy_temp c1998.std_gy_temp c1999.std_gy_temp c2000.std_gy_temp c2001.std_gy_temp ...
    c2002.std_gy_temp c2003.std_gy_temp c2004.std_gy_temp c2005.std_gy_temp c2006.std_gy_temp c2007.std_gy_temp ...
    c2008.std_gy_temp c2009.std_gy_temp c2010.std_gy_temp];
mer_std_salt = [c1997.std_gy_salt c1998.std_gy_salt c1999.std_gy_salt c2000.std_gy_salt c2001.std_gy_salt ...
    c2002.std_gy_salt c2003.std_gy_salt c2004.std_gy_salt c2005.std_gy_salt c2006.std_gy_salt c2007.std_gy_salt ...
    c2008.std_gy_salt c2009.std_gy_salt c2010.std_gy_salt];
%bot merged
mer_std_no3_b = [c1997.std_gy_no3_b c1998.std_gy_no3_b c1999.std_gy_no3_b c2000.std_gy_no3_b c2001.std_gy_no3_b ...
    c2002.std_gy_no3_b c2003.std_gy_no3_b c2004.std_gy_no3_b c2005.std_gy_no3_b c2006.std_gy_no3_b c2007.std_gy_no3_b ...
    c2008.std_gy_no3_b c2009.std_gy_no3_b c2010.std_gy_no3_b];
mer_std_nh4_b = [c1997.std_gy_nh4_b c1998.std_gy_nh4_b c1999.std_gy_nh4_b c2000.std_gy_nh4_b c2001.std_gy_nh4_b ...
    c2002.std_gy_nh4_b c2003.std_gy_nh4_b c2004.std_gy_nh4_b c2005.std_gy_nh4_b c2006.std_gy_nh4_b c2007.std_gy_nh4_b ...
    c2008.std_gy_nh4_b c2009.std_gy_nh4_b c2010.std_gy_nh4_b];
mer_std_chl_b = [c1997.std_gy_chl_b c1998.std_gy_chl_b c1999.std_gy_chl_b c2000.std_gy_chl_b c2001.std_gy_chl_b ...
    c2002.std_gy_chl_b c2003.std_gy_chl_b c2004.std_gy_chl_b c2005.std_gy_chl_b c2006.std_gy_chl_b c2007.std_gy_chl_b ...
    c2008.std_gy_chl_b c2009.std_gy_chl_b c2010.std_gy_chl_b];
mer_std_temp_b = [c1997.std_gy_temp_b c1998.std_gy_temp_b c1999.std_gy_temp_b c2000.std_gy_temp_b c2001.std_gy_temp_b ...
    c2002.std_gy_temp_b c2003.std_gy_temp_b c2004.std_gy_temp_b c2005.std_gy_temp_b c2006.std_gy_temp_b c2007.std_gy_temp_b ...
    c2008.std_gy_temp_b c2009.std_gy_temp_b c2010.std_gy_temp_b];
mer_std_salt_b = [c1997.std_gy_salt_b c1998.std_gy_salt_b c1999.std_gy_salt_b c2000.std_gy_salt_b c2001.std_gy_salt_b ...
    c2002.std_gy_salt_b c2003.std_gy_salt_b c2004.std_gy_salt_b c2005.std_gy_salt_b c2006.std_gy_salt_b c2007.std_gy_salt_b ...
    c2008.std_gy_salt_b c2009.std_gy_salt_b c2010.std_gy_salt_b];

com_eom = [1:length(1997:2010)*12];

% % for i = 1:length(t_tIc)+length(t_tIc2)
% %     if i < 13
% %         com_t_tic{i}=t_tIc{i}
% %     else
% %         com_t_tic{i}=t_tIc2{i-12}
% %     end
% % end

% 1997~2018 : 22yr
% 1st regime : 1997 ~ 2003.12 (from koem nh4 shift)
% 2nd regime : 2004.01 ~ 2009.12 (from koem no3 shift)
% 3rd regime : 2010.01 ~ 2018.12 (from koem no3 shift)

% 1st
clim1.ref_date(84,:) % 1st regime : 1997 ~ 2003.12 (from koem nh4 shift)
% 2nd
clim1.ref_date(85,:) % 2nd regime : 2004.01 ~ 2009.12 (from koem no3 shift)
clim1.ref_date(156,:)
% 3rd
clim1.ref_date(157,:) % 3rd regime : 2010.01 ~ 2018.12 (from koem no3 shift)
clim1.ref_date(end,:)

return
for i =  1:12
%% model
    % 1st
    regm_1st_no3(i)=mean(mer_no3(i:12:84));
    regm_1st_nh4(i)=mean(mer_nh4(i:12:84));
    regm_1st_chl(i)=mean(mer_chl(i:12:84));
    regm_1st_temp(i)=mean(mer_temp(i:12:84));
    regm_1st_salt(i)=mean(mer_salt(i:12:84));
    % 2nd
    regm_2nd_no3(i)=mean(mer_no3(84+i:12:156));
    regm_2nd_nh4(i)=mean(mer_nh4(84+i:12:156));
    regm_2nd_chl(i)=mean(mer_chl(84+i:12:156));
    regm_2nd_temp(i)=mean(mer_temp(84+i:12:156));
    regm_2nd_salt(i)=mean(mer_salt(84+i:12:156));  
  
    % 1st
    regm_1st_no3_b(i)=mean(mer_no3_b(i:12:84));
    regm_1st_nh4_b(i)=mean(mer_nh4_b(i:12:84));
    regm_1st_chl_b(i)=mean(mer_chl_b(i:12:84));
    regm_1st_temp_b(i)=mean(mer_temp_b(i:12:84));
    regm_1st_salt_b(i)=mean(mer_salt_b(i:12:84));
    % 2nd
    regm_2nd_no3_b(i)=mean(mer_no3_b(84+i:12:156));
    regm_2nd_nh4_b(i)=mean(mer_nh4_b(84+i:12:156));
    regm_2nd_chl_b(i)=mean(mer_chl_b(84+i:12:156));
    regm_2nd_temp_b(i)=mean(mer_temp_b(84+i:12:156));
    regm_2nd_salt_b(i)=mean(mer_salt_b(84+i:12:156)); 
%std
    % 1st
    clearvars size2
    size2=size(mer_sp_no3(:,i:12:84),2);
    regm_std_1st_no3(i)=std(reshape(mer_sp_no3(:,i:12:84),1,9*size2));
    regm_std_1st_nh4(i)=std(reshape(mer_sp_nh4(:,i:12:84),1,9*size2));
    regm_std_1st_chl(i)=std(reshape(mer_sp_chl(:,i:12:84),1,9*size2));
    regm_std_1st_temp(i)=std(reshape(mer_sp_temp(:,i:12:84),1,9*size2));
    regm_std_1st_salt(i)=std(reshape(mer_sp_salt(:,i:12:84),1,9*size2));
    % 2nd
    clearvars size2
    size2=size(mer_sp_no3(:,84+i:12:156),2);
    regm_std_2nd_no3(i)=std(reshape(mer_sp_no3(:,84+i:12:156),1,9*size2));
    regm_std_2nd_nh4(i)=std(reshape(mer_sp_nh4(:,84+i:12:156),1,9*size2));
    regm_std_2nd_chl(i)=std(reshape(mer_sp_chl(:,84+i:12:156),1,9*size2));
    regm_std_2nd_temp(i)=std(reshape(mer_sp_temp(:,84+i:12:156),1,9*size2));
    regm_std_2nd_salt(i)=std(reshape(mer_sp_salt(:,84+i:12:156),1,9*size2));  
  
    % 1st
    clearvars size2
    size2=size(mer_sp_no3(:,i:12:84),2);
    regm_std_1st_no3_b(i)=std(reshape(mer_sp_no3_b(:,i:12:84),1,9*size2));
    regm_std_1st_nh4_b(i)=std(reshape(mer_sp_nh4_b(:,i:12:84),1,9*size2));
    regm_std_1st_chl_b(i)=std(reshape(mer_sp_chl_b(:,i:12:84),1,9*size2));
    regm_std_1st_temp_b(i)=std(reshape(mer_sp_temp_b(:,i:12:84),1,9*size2));
    regm_std_1st_salt_b(i)=std(reshape(mer_sp_salt_b(:,i:12:84),1,9*size2));
    % 2nd
    clearvars size2
    size2=size(mer_sp_no3(:,84+i:12:156),2);
    regm_std_2nd_no3_b(i)=std(reshape(mer_sp_no3_b(:,84+i:12:156),1,9*size2));
    regm_std_2nd_nh4_b(i)=std(reshape(mer_sp_nh4_b(:,84+i:12:156),1,9*size2));
    regm_std_2nd_chl_b(i)=std(reshape(mer_sp_chl_b(:,84+i:12:156),1,9*size2));
    regm_std_2nd_temp_b(i)=std(reshape(mer_sp_temp_b(:,84+i:12:156),1,9*size2));
    regm_std_2nd_salt_b(i)=std(reshape(mer_sp_salt_b(:,84+i:12:156),1,9*size2)); 
    
    
%% obs
    % 1st
    obs_regm_1st_no3(i)=nanmean(clim1.obm_gy_no3(i:12:84));
    obs_regm_1st_nh4(i)=nanmean(clim1.obm_gy_nh4(i:12:84));
    obs_regm_1st_chl(i)=nanmean(clim1.obm_gy_chl(i:12:84));
    obs_regm_1st_temp(i)=nanmean(clim1.obm_gy_temp(i:12:84));
    obs_regm_1st_salt(i)=nanmean(clim1.obm_gy_salt(i:12:84));
    % 2nd
    obs_regm_2nd_no3(i)=nanmean(clim1.obm_gy_no3(84+i:12:156));
    obs_regm_2nd_nh4(i)=nanmean(clim1.obm_gy_nh4(84+i:12:156));
    obs_regm_2nd_chl(i)=nanmean(clim1.obm_gy_chl(84+i:12:156));
    obs_regm_2nd_temp(i)=nanmean(clim1.obm_gy_temp(84+i:12:156));
    obs_regm_2nd_salt(i)=nanmean(clim1.obm_gy_salt(84+i:12:156));  
  
    % 1st
    obs_regm_1st_no3_b(i)=nanmean(clim1.obm_gy_no3_b(i:12:84));
    obs_regm_1st_nh4_b(i)=nanmean(clim1.obm_gy_nh4_b(i:12:84));
    obs_regm_1st_chl_b(i)=nanmean(clim1.obm_gy_chl_b(i:12:84));
    obs_regm_1st_temp_b(i)=nanmean(clim1.obm_gy_temp_b(i:12:84));
    obs_regm_1st_salt_b(i)=nanmean(clim1.obm_gy_salt_b(i:12:84));
    % 2nd
    obs_regm_2nd_no3_b(i)=nanmean(clim1.obm_gy_no3_b(84+i:12:156));
    obs_regm_2nd_nh4_b(i)=nanmean(clim1.obm_gy_nh4_b(84+i:12:156));
    obs_regm_2nd_chl_b(i)=nanmean(clim1.obm_gy_chl_b(84+i:12:156));
    obs_regm_2nd_temp_b(i)=nanmean(clim1.obm_gy_temp_b(84+i:12:156));
    obs_regm_2nd_salt_b(i)=nanmean(clim1.obm_gy_salt_b(84+i:12:156)); 
% std    
    % 1st
    obs_std_regm_1st_no3(i)=nanstd(clim1.obm_gy_no3(i:12:84));
    obs_std_regm_1st_nh4(i)=nanstd(clim1.obm_gy_nh4(i:12:84));
    obs_std_regm_1st_chl(i)=nanstd(clim1.obm_gy_chl(i:12:84));
    obs_std_regm_1st_temp(i)=nanstd(clim1.obm_gy_temp(i:12:84));
    obs_std_regm_1st_salt(i)=nanstd(clim1.obm_gy_salt(i:12:84));
    % 2nd
    obs_std_regm_2nd_no3(i)=nanstd(clim1.obm_gy_no3(84+i:12:156));
    obs_std_regm_2nd_nh4(i)=nanstd(clim1.obm_gy_nh4(84+i:12:156));
    obs_std_regm_2nd_chl(i)=nanstd(clim1.obm_gy_chl(84+i:12:156));
    obs_std_regm_2nd_temp(i)=nanstd(clim1.obm_gy_temp(84+i:12:156));
    obs_std_regm_2nd_salt(i)=nanstd(clim1.obm_gy_salt(84+i:12:156));  
  
    % 1st
    obs_std_regm_1st_no3_b(i)=nanstd(clim1.obm_gy_no3_b(i:12:84));
    obs_std_regm_1st_nh4_b(i)=nanstd(clim1.obm_gy_nh4_b(i:12:84));
    obs_std_regm_1st_chl_b(i)=nanstd(clim1.obm_gy_chl_b(i:12:84));
    obs_std_regm_1st_temp_b(i)=nanstd(clim1.obm_gy_temp_b(i:12:84));
    obs_std_regm_1st_salt_b(i)=nanstd(clim1.obm_gy_salt_b(i:12:84));
    % 2nd
    obs_std_regm_2nd_no3_b(i)=nanstd(clim1.obm_gy_no3_b(84+i:12:156));
    obs_std_regm_2nd_nh4_b(i)=nanstd(clim1.obm_gy_nh4_b(84+i:12:156));
    obs_std_regm_2nd_chl_b(i)=nanstd(clim1.obm_gy_chl_b(84+i:12:156));
    obs_std_regm_2nd_temp_b(i)=nanstd(clim1.obm_gy_temp_b(84+i:12:156));
    obs_std_regm_2nd_salt_b(i)=nanstd(clim1.obm_gy_salt_b(84+i:12:156));      
end


clearvars com_t_tic
% com_t_tic = {'2001.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             'c2002.01','02','03','04','05','06','07','08','09','10','11','12', ...
%             '2003.01','02','03','04','05','06','07','08','09','10','11','12',...}
%             '2004.01','02','03','04','05','06','07','08','09','10','11','12'}

        
com_t_tic = {'1997','1998','1999','2000','2001', ...
            '2002', ...
            '2003',...
            '2004',...
            '2005',...
            '2006',...
            '2007',...
            '2008', '2009','2010'}
        
           
%% GY plot

% obs_tdx =12*(year_f - 1997)+1:12*(year_f - 1997)+12;

% 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_no3 + regm_std_1st_no3;
        lower_bound_plt = regm_1st_no3 - regm_std_1st_no3;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_no3./14,obs_std_regm_2nd_no3./14,'color','r','linew',1.5);
        plot(obs_regm_1st_no3./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL NO3']);
        xlabel('time(month)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 70])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_no3'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_no3,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_no3)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_no3 + regm_std_2nd_no3;
        lower_bound_plt = regm_2nd_no3 - regm_std_2nd_no3;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_no3./14,obs_std_regm_2nd_no3./14,'color','r','linew',1.5);
        plot(obs_regm_2nd_no3./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL NO3']);
        xlabel('time(month)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 70])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_no3'),'-dpng')
        
%         

% 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_nh4 + regm_std_1st_nh4;
        lower_bound_plt = regm_1st_nh4 - regm_std_1st_nh4;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_nh4./14,obs_std_regm_2nd_nh4./14,'color','r','linew',1.5);
        plot(obs_regm_1st_nh4./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_nh4'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_nh4,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_nh4)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_nh4 + regm_std_2nd_nh4;
        lower_bound_plt = regm_2nd_nh4 - regm_std_2nd_nh4;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_nh4./14,obs_std_regm_2nd_nh4./14,'color','r','linew',1.5);
        plot(obs_regm_2nd_nh4./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL nh4']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_nh4'),'-dpng')
        
        
% 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_chl + regm_std_1st_chl;
        lower_bound_plt = regm_1st_chl - regm_std_1st_chl;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_chl,obs_std_regm_2nd_chl,'color','r','linew',1.5);
        plot(obs_regm_1st_chl,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_chl'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_chl,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_chl)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_chl + regm_std_2nd_chl;
        lower_bound_plt = regm_2nd_chl - regm_std_2nd_chl;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_chl,obs_std_regm_2nd_chl,'color','r','linew',1.5);
        plot(obs_regm_2nd_chl,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_chl'),'-dpng')
        
        
 % 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_temp + regm_std_1st_temp;
        lower_bound_plt = regm_1st_temp - regm_std_1st_temp;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_temp,obs_std_regm_2nd_temp,'color','r','linew',1.5);
        plot(obs_regm_1st_temp,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL temp']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_temp'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_temp,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_temp)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_temp + regm_std_2nd_temp;
        lower_bound_plt = regm_2nd_temp - regm_std_2nd_temp;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_temp,obs_std_regm_2nd_temp,'color','r','linew',1.5);
        plot(obs_regm_2nd_temp,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL temp']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_temp'),'-dpng')
        
        
% 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_salt + regm_std_1st_salt;
        lower_bound_plt = regm_1st_salt - regm_std_1st_salt;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_salt,obs_std_regm_2nd_salt,'color','r','linew',1.5);
        plot(obs_regm_1st_salt,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL salt']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_salt'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_salt,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_salt)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_salt + regm_std_2nd_salt;
        lower_bound_plt = regm_2nd_salt - regm_std_2nd_salt;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_salt,obs_std_regm_2nd_salt,'color','r','linew',1.5);
        plot(obs_regm_2nd_salt,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL salt']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_salt'),'-dpng')
        
        %% bot
 
 % 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_no3_b + regm_std_1st_no3_b;
        lower_bound_plt = regm_1st_no3_b - regm_std_1st_no3_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_no3_b./14,obs_std_regm_2nd_no3_b./14,'color','r','linew',1.5);
        plot(obs_regm_1st_no3_b./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL no3 bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 70])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_no3_b'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_no3_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_no3_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_no3_b + regm_std_2nd_no3_b;
        lower_bound_plt = regm_2nd_no3_b - regm_std_2nd_no3_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_no3_b./14,obs_std_regm_2nd_no3_b./14,'color','r','linew',1.5);
        plot(obs_regm_2nd_no3_b./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL no3 bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('no3 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 70])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_no3_b'),'-dpng')
        
%         

% 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_nh4_b + regm_std_1st_nh4_b;
        lower_bound_plt = regm_1st_nh4_b - regm_std_1st_nh4_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_nh4_b./14,obs_std_regm_2nd_nh4_b./14,'color','r','linew',1.5);
        plot(obs_regm_1st_nh4_b./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL nh4 bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_nh4_b'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_nh4_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_nh4_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_nh4_b + regm_std_2nd_nh4_b;
        lower_bound_plt = regm_2nd_nh4_b - regm_std_2nd_nh4_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_nh4_b./14,obs_std_regm_2nd_nh4_b./14,'color','r','linew',1.5);
        plot(obs_regm_2nd_nh4_b./14,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL nh4 bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('nh4 (umol/m^3)','fontsize',13)
        grid on
        ylim([0 14])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_nh4_b'),'-dpng')
        
        
% 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_chl_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_chl_b + regm_std_1st_chl_b;
        lower_bound_plt = regm_1st_chl_b - regm_std_1st_chl_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_chl_b,obs_std_regm_2nd_chl_b,'color','r','linew',1.5);
        plot(obs_regm_1st_chl_b,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL chl bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_chl_b'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_chl_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_chl_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_chl_b + regm_std_2nd_chl_b;
        lower_bound_plt = regm_2nd_chl_b - regm_std_2nd_chl_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_chl_b,obs_std_regm_2nd_chl_b,'color','r','linew',1.5);
        plot(obs_regm_2nd_chl_b,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL chl bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        grid on
        ylim([0 15])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_chl_b'),'-dpng')
        
        
 % 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_temp_b + regm_std_1st_temp_b;
        lower_bound_plt = regm_1st_temp_b - regm_std_1st_temp_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_temp_b,obs_std_regm_2nd_temp_b,'color','r','linew',1.5);
        plot(obs_regm_1st_temp_b,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL temp bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_temp_b'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_temp_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_temp_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_temp_b + regm_std_2nd_temp_b;
        lower_bound_plt = regm_2nd_temp_b - regm_std_2nd_temp_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_temp_b,obs_std_regm_2nd_temp_b,'color','r','linew',1.5);
        plot(obs_regm_2nd_temp_b,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL temp bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_temp_b'),'-dpng')
        
        
% 1st
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_1st_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_1st_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_1st_salt_b + regm_std_1st_salt_b;
        lower_bound_plt = regm_1st_salt_b - regm_std_1st_salt_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_1st_salt_b,obs_std_regm_2nd_salt_b,'color','r','linew',1.5);
        plot(obs_regm_1st_salt_b,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 1997~2003 monthly KOEM OBS vs. MODEL salt bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2003_month_KOEM_OBS_vs_MODEL_salt_b'),'-dpng')
 % 2nd
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(regm_2nd_salt_b,'b','linew',2);
        clearvars *_bound_plt nonan_data_plt discon_p
        nonan_data_plt = find(isnan(regm_2nd_salt_b)==0);
        nonan_diff=nonan_data_plt(2:end) - nonan_data_plt(1:end-1);
        discon_p = find(nonan_diff ~= 1); % find discontinuous point
        discon_p(end+1) = length(nonan_data_plt); %end point_ %end point
        upper_bound_plt = regm_2nd_salt_b + regm_std_2nd_salt_b;
        lower_bound_plt = regm_2nd_salt_b - regm_std_2nd_salt_b;
        
%         for i = 1:length(upper_bound_plt)
%         upper_bnd(i)=max(upper_bound_plt(i),lower_bound_plt(i));
%         lower_bnd(i)=min(upper_bound_plt(i),lower_bound_plt(i));
%         end
        
        for i = 1:length(discon_p)
        if i == 1
        patch([nonan_data_plt(1:discon_p(1)) fliplr(nonan_data_plt(1:discon_p(1)))], [upper_bound_plt(nonan_data_plt(1:discon_p(1))) fliplr(lower_bound_plt(nonan_data_plt(1:discon_p(1))))],[0.5 0.5 0.5]);
        else
        patch([nonan_data_plt(discon_p(i-1)+1:discon_p(i)) fliplr(nonan_data_plt(discon_p(i-1)+1:discon_p(i)))], [upper_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))) fliplr(lower_bound_plt(nonan_data_plt(discon_p(i-1)+1:discon_p(i))))],'b');
        end
        end 
        errorbar(obs_regm_2nd_salt_b,obs_std_regm_2nd_salt_b,'color','r','linew',1.5);
        plot(obs_regm_2nd_salt_b,'r+','linew',1.5); xlim([1 12]);
        alpha(0.3) %transparency
        
        title(['GY 2004~2009 monthly KOEM OBS vs. MODEL salt bot.']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        ylim([0 35])
        set(gca,'fontsize',13)
        print(fig,strcat('2004~2009_month_KOEM_OBS_vs_MODEL_salt_b'),'-dpng')
        
return

%% GY plot
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_po4,'b','linew',2);
%         plot(obm_gy_po4./14,'r-','linew',2); 
        xlim([1 12]);
%         plot(1:length(obm_gy_po4),obm_gy_po4./14 + obs_std_gy_po4./14,'m-','linew',2);
%         plot(1:length(obm_gy_po4),obm_gy_po4./14 - obs_std_gy_po4./14,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL po4']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_po4_gy'),'-dpng')


fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3,'b','linew',2);
        plot(obm_gy_no3./14,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_no3),obm_gy_no3./14 + obs_std_gy_no3./14,'m-','linew',2);
        plot(1:length(obm_gy_no3),obm_gy_no3./14 - obs_std_gy_no3./14,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL NO3']);
        xlabel('time(month)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        ylim([0 55]);
        grid on;
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_no3_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4,'b','linew',2);
        plot(obm_gy_nh4./14,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_nh4),obm_gy_nh4./14 + obs_std_gy_nh4./14,'m-','linew',2);
        plot(1:length(obm_gy_nh4),obm_gy_nh4./14 - obs_std_gy_nh4./14,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL NH4']);
        xlabel('time(month)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        ylim([0 14]);
        grid on;
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_NH4_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl,'b','linew',2);
        plot(obm_gy_chl,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_chl),obm_gy_chl + obs_std_gy_chl,'m-','linew',2);
        plot(1:length(obm_gy_chl),obm_gy_chl - obs_std_gy_chl,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL chl']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        ylim([0 15])
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_chl_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp,'b','linew',2);
        plot(obm_gy_temp,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_temp),obm_gy_temp + obs_std_gy_temp,'m-','linew',2);
        plot(1:length(obm_gy_temp),obm_gy_temp - obs_std_gy_temp,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL temp']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_temp_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt,'b','linew',2);
        plot(obm_gy_salt,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_salt),obm_gy_salt + obs_std_gy_salt,'m-','linew',2);
        plot(1:length(obm_gy_salt),obm_gy_salt - obs_std_gy_salt,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL salt']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_salt_gy'),'-dpng')            
 %% bot
  fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_po4_b,'b','linew',2);
%         plot(obm_gy_po4_b./14,'r-','linew',2);
        xlim([1 12]);
%         plot(1:length(obm_gy_po4_b),obm_gy_po4_b./14 + obs_std_gy_po4_b./14,'m-','linew',2);
%         plot(1:length(obm_gy_po4_b),obm_gy_po4_b./14 - obs_std_gy_po4_b./14,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL po4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('po4 (umol/m^3)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_po4_bot_gy'),'-dpng')
 
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_no3_b,'b','linew',2);
        plot(obm_gy_no3_b./14,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 + obs_std_gy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_gy_no3_b),obm_gy_no3_b./14 - obs_std_gy_no3_b./14,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('NO3 (umol/m^3)','fontsize',13)
        ylim([0 55]);
        grid on;
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_no3_bot_gy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_nh4_b,'b','linew',2);
        plot(obm_gy_nh4_b./14,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 + obs_std_gy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_gy_nh4_b),obm_gy_nh4_b./14 - obs_std_gy_nh4_b./14,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('NH4 (umol/m^3)','fontsize',13)
        ylim([0 14])
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_NH4_bot_gy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_chl_b,'b','linew',2);
        plot(obm_gy_chl_b,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_chl_b),obm_gy_chl_b + obs_std_gy_chl_b,'m-','linew',2);
        plot(1:length(obm_gy_chl_b),obm_gy_chl_b - obs_std_gy_chl_b,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL chl bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('chl (ug/L)','fontsize',13)
        ylim([0 15])
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_chl_bot_gy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_temp_b,'b','linew',2);
        plot(obm_gy_temp_b,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_temp_b),obm_gy_temp_b + obs_std_gy_temp_b,'m-','linew',2);
        plot(1:length(obm_gy_temp_b),obm_gy_temp_b - obs_std_gy_temp_b,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL temp bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('temp (o^C)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_temp_bot_gy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g--','linew',2); % zero lines
        plot(gy_salt_b,'b','linew',2);
        plot(obm_gy_salt_b,'r-','linew',2); xlim([1 12]);
        plot(1:length(obm_gy_salt_b),obm_gy_salt_b + obs_std_gy_salt_b,'m-','linew',2);
        plot(1:length(obm_gy_salt_b),obm_gy_salt_b - obs_std_gy_salt_b,'m-','linew',2);
        title(['GY 1997~2010 monthly KOEM OBS vs. MODEL salt bot']);
        xlabel('time(month)','fontsize',13)
        ylabel('salt (psu)','fontsize',13)
        grid on
        set(gca,'fontsize',13)
        print(fig,strcat('1997~2010_month_KOEM_OBS_vs_MODEL_salt_bot_gy'),'-dpng') 
        
return        
        
%%south gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_no3,'b','linew',2);
        plot(pcase.sgy_no3,'k','linew',2);
        plot(mp_p_ws.sgy_no3,'g','linew',2); 
        plot(obm_sgy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 + obs_std_sgy_no3./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3),obm_sgy_no3./14 - obs_std_sgy_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_no3_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_nh4,'b','linew',2);
        plot(pcase.sgy_nh4,'k','linew',2);
        plot(mp_p_ws.sgy_nh4,'g','linew',2);
        plot(obm_sgy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 + obs_std_sgy_nh4./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4),obm_sgy_nh4./14 - obs_std_sgy_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_NH4_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_chl,'b','linew',2);
        plot(pcase.sgy_chl,'k','linew',2);
        plot(mp_p_ws.sgy_chl,'g','linew',2);
        plot(obm_sgy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl),obm_sgy_chl + obs_std_sgy_chl,'m-','linew',2);
        plot(1:length(obm_sgy_chl),obm_sgy_chl - obs_std_sgy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_chl_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_temp,'b','linew',2);
        plot(pcase.sgy_temp,'k','linew',2);
        plot(mp_p_ws.sgy_temp,'g','linew',2);
        plot(obm_sgy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp),obm_sgy_temp + obs_std_sgy_temp,'m-','linew',2);
        plot(1:length(obm_sgy_temp),obm_sgy_temp - obs_std_sgy_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_temp_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt,'b','linew',2);
        plot(pcase.sgy_salt,'k','linew',2);
        plot(mp_p_ws.sgy_salt,'g','linew',2);   
        plot(obm_sgy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt),obm_sgy_salt + obs_std_sgy_salt,'m-','linew',2);
        plot(1:length(obm_sgy_salt),obm_sgy_salt - obs_std_sgy_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_salt_south_sgy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_no3_b,'b','linew',2);
        plot(pcase.sgy_no3_b,'k','linew',2);
        plot(mp_p_ws.sgy_no3_b,'g','linew',2);
        plot(obm_sgy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 + obs_std_sgy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_no3_b),obm_sgy_no3_b./14 - obs_std_sgy_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_no3_bot_south_sgy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_nh4_b,'b','linew',2);
        plot(pcase.sgy_nh4_b,'k','linew',2);
        plot(mp_p_ws.sgy_nh4_b,'g','linew',2); 
        plot(obm_sgy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 + obs_std_sgy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_sgy_nh4_b),obm_sgy_nh4_b./14 - obs_std_sgy_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_NH4_bot_south_sgy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_chl_b,'b','linew',2);
        plot(pcase.sgy_chl_b,'k','linew',2);
        plot(mp_p_ws.sgy_chl_b,'g','linew',2);  
        plot(obm_sgy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b + obs_std_sgy_chl_b,'m-','linew',2);
        plot(1:length(obm_sgy_chl_b),obm_sgy_chl_b - obs_std_sgy_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_chl_bot_south_sgy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_temp_b,'b','linew',2);
        plot(pcase.sgy_temp_b,'k','linew',2);
        plot(mp_p_ws.sgy_temp_b,'g','linew',2);        
        plot(obm_sgy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b + obs_std_sgy_temp_b,'m-','linew',2);
        plot(1:length(obm_sgy_temp_b),obm_sgy_temp_b - obs_std_sgy_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_temp_bot_south_sgy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.sgy_salt_b,'b','linew',2);
        plot(pcase.sgy_salt_b,'k','linew',2);
        plot(mp_p_ws.sgy_salt_b,'g','linew',2);    
        plot(obm_sgy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b + obs_std_sgy_salt_b,'m-','linew',2);
        plot(1:length(obm_sgy_salt_b),obm_sgy_salt_b - obs_std_sgy_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_salt_bot_south_sgy'),'-dpng')    
        
%%east gy plot

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_no3,'b','linew',2);
        plot(pcase.egy_no3,'k','linew',2);
        plot(mp_p_ws.egy_no3,'g','linew',2);  
        plot(obm_egy_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 + obs_std_egy_no3./14,'m-','linew',2);
        plot(1:length(obm_egy_no3),obm_egy_no3./14 - obs_std_egy_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_no3_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines       
        plot(mpdata.egy_nh4,'b','linew',2);
        plot(pcase.egy_nh4,'k','linew',2);
        plot(mp_p_ws.egy_nh4,'g','linew',2); 
        plot(obm_egy_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 + obs_std_egy_nh4./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4),obm_egy_nh4./14 - obs_std_egy_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_NH4_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_chl,'b','linew',2);
        plot(pcase.egy_chl,'k','linew',2);
        plot(mp_p_ws.egy_chl,'g','linew',2); 
        plot(obm_egy_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl),obm_egy_chl + obs_std_egy_chl,'m-','linew',2);
        plot(1:length(obm_egy_chl),obm_egy_chl - obs_std_egy_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_chl_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp,'b','linew',2);
        plot(pcase.egy_temp,'k','linew',2);
        plot(mp_p_ws.egy_temp,'g','linew',2); 
        plot(obm_egy_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp),obm_egy_temp + obs_std_egy_temp,'m-','linew',2);
        plot(1:length(obm_egy_temp),obm_egy_temp - obs_std_egy_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_temp_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt,'b','linew',2);
        plot(pcase.egy_salt,'k','linew',2);
        plot(mp_p_ws.egy_salt,'g','linew',2); 
        plot(obm_egy_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt),obm_egy_salt + obs_std_egy_salt,'m-','linew',2);
        plot(1:length(obm_egy_salt),obm_egy_salt - obs_std_egy_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_salt_egy'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_no3_b,'b','linew',2);
        plot(pcase.egy_no3_b,'k','linew',2);
        plot(mp_p_ws.egy_no3_b,'g','linew',2); 
        plot(obm_egy_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 + obs_std_egy_no3_b./14,'m-','linew',2);
        plot(1:length(obm_egy_no3_b),obm_egy_no3_b./14 - obs_std_egy_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_no3_bot_egy'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_nh4_b,'b','linew',2);
        plot(pcase.egy_nh4_b,'k','linew',2);
        plot(mp_p_ws.egy_nh4_b,'g','linew',2); 
        plot(obm_egy_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 + obs_std_egy_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_egy_nh4_b),obm_egy_nh4_b./14 - obs_std_egy_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_NH4_bot_egy'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_chl_b,'b','linew',2);
        plot(pcase.egy_chl_b,'k','linew',2);
        plot(mp_p_ws.egy_chl_b,'g','linew',2); 
        plot(obm_egy_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b + obs_std_egy_chl_b,'m-','linew',2);
        plot(1:length(obm_egy_chl_b),obm_egy_chl_b - obs_std_egy_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_chl_bot_egy'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_temp_b,'b','linew',2);
        plot(pcase.egy_temp_b,'k','linew',2);
        plot(mp_p_ws.egy_temp_b,'g','linew',2); 
        plot(obm_egy_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b + obs_std_egy_temp_b,'m-','linew',2);
        plot(1:length(obm_egy_temp_b),obm_egy_temp_b - obs_std_egy_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_temp_bot_egy'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.egy_salt_b,'b','linew',2);
        plot(pcase.egy_salt_b,'k','linew',2);
        plot(mp_p_ws.egy_salt_b,'g','linew',2); 
        plot(obm_egy_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b + obs_std_egy_salt_b,'m-','linew',2);
        plot(1:length(obm_egy_salt_b),obm_egy_salt_b - obs_std_egy_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_salt_bot_egy'),'-dpng')    
        
%% jinju

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_no3,'b','linew',2);
        plot(pcase.jj_no3,'k','linew',2);
        plot(mp_p_ws.jj_no3,'g','linew',2); 
        plot(obm_jj_no3./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 + obs_std_jj_no3./14,'m-','linew',2);
        plot(1:length(obm_jj_no3),obm_jj_no3./14 - obs_std_jj_no3./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_no3_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4,'b','linew',2);
        plot(pcase.jj_nh4,'k','linew',2);
        plot(mp_p_ws.jj_nh4,'g','linew',2); 
        plot(obm_jj_nh4./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 + obs_std_jj_nh4./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4),obm_jj_nh4./14 - obs_std_jj_nh4./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_NH4_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_chl,'b','linew',2);
        plot(pcase.jj_chl,'k','linew',2);
        plot(mp_p_ws.jj_chl,'g','linew',2); 
        plot(obm_jj_chl,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl),obm_jj_chl + obs_std_jj_chl,'m-','linew',2);
        plot(1:length(obm_jj_chl),obm_jj_chl - obs_std_jj_chl,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_chl_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp,'b','linew',2);
        plot(pcase.jj_temp,'k','linew',2);
        plot(mp_p_ws.jj_temp,'g','linew',2); 
        plot(obm_jj_temp,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp),obm_jj_temp + obs_std_jj_temp,'m-','linew',2);
        plot(1:length(obm_jj_temp),obm_jj_temp - obs_std_jj_temp,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_temp_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt,'b','linew',2);
        plot(pcase.jj_salt,'k','linew',2);
        plot(mp_p_ws.jj_salt,'g','linew',2); 
        plot(obm_jj_salt,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt),obm_jj_salt + obs_std_jj_salt,'m-','linew',2);
        plot(1:length(obm_jj_salt),obm_jj_salt - obs_std_jj_salt,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_salt_jj'),'-dpng')            
 %% bot
 
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_no3_b,'b','linew',2);
        plot(pcase.jj_no3_b,'k','linew',2);
        plot(mp_p_ws.jj_no3_b,'g','linew',2); 
        plot(obm_jj_no3_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 + obs_std_jj_no3_b./14,'m-','linew',2);
        plot(1:length(obm_jj_no3_b),obm_jj_no3_b./14 - obs_std_jj_no3_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NO3 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NO3 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_no3_bot_jj'),'-dpng')

fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_nh4_b,'b','linew',2);
        plot(pcase.jj_nh4_b,'k','linew',2);
        plot(mp_p_ws.jj_nh4_b,'g','linew',2); 
        plot(obm_jj_nh4_b./14,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 + obs_std_jj_nh4_b./14,'m-','linew',2);
        plot(1:length(obm_jj_nh4_b),obm_jj_nh4_b./14 - obs_std_jj_nh4_b./14,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL NH4 bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('NH4 (umol/m^3)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_NH4_bot_jj'),'-dpng')
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_chl_b,'b','linew',2);
        plot(pcase.jj_chl_b,'k','linew',2);
        plot(mp_p_ws.jj_chl_b,'g','linew',2); 
        plot(obm_jj_chl_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b + obs_std_jj_chl_b,'m-','linew',2);
        plot(1:length(obm_jj_chl_b),obm_jj_chl_b - obs_std_jj_chl_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL chl bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('chl (ug/L)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold'); ylim([0 inf])
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_chl_bot_jj'),'-dpng') 
        
        
 fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_temp_b,'b','linew',2);
        plot(pcase.jj_temp_b,'k','linew',2);
        plot(mp_p_ws.jj_temp_b,'g','linew',2); 
        plot(obm_jj_temp_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b + obs_std_jj_temp_b,'m-','linew',2);
        plot(1:length(obm_jj_temp_b),obm_jj_temp_b - obs_std_jj_temp_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL temp bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('temp (o^C)','fontsize',16)
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_temp_bot_jj'),'-dpng')    
        
        
fig = figure; hold on;
%         plot(zeros(365,1),'g','linew',2); % zero lines
        plot(mpdata.jj_salt_b,'b','linew',2);
        plot(pcase.jj_salt_b,'k','linew',2);
        plot(mp_p_ws.jj_salt_b,'g','linew',2); 
        plot(obm_jj_salt_b,'r-','linew',2); xlim([1 365]);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b + obs_std_jj_salt_b,'m-','linew',2);
        plot(1:length(obm_jj_salt_b),obm_jj_salt_b - obs_std_jj_salt_b,'m-','linew',2);
        title(['GY 2001 daily KOEM OBS vs. MODEL salt bot']);
        xlabel('time(days on 2001)','fontsize',16)
        ylabel('salt (psu)','fontsize',16)
        xticks([0:50:366])
        grid on
        set(gca,'fontsize',16,'fontweight','bold')
        print(fig,strcat('1997~2010_daily_KOEM_OBS_vs_MODEL_salt_bot_jj'),'-dpng')    


