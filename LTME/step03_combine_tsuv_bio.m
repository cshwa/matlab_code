clc;clear all;close all
% 
% for i=2011:2017
%     eval(['!copy roms_river_auto_new8_',num2str(i),'_po4_8rivers.nc roms_river_auto_new8_',num2str(i),'_yangtze80s_po4_8rivers.nc'])
% end

year=2011;

% tsuv_file=['roms_river_auto_new8_',num2str(year),'_po4_8rivers.nc'];
% bio_file=['roms_river_auto_new8_',num2str(year),'_fennel_po4_8rivers.nc'];

% tsuv_file=['roms_river_auto_new8_',num2str(year),'_yangtze_po4_8rivers.nc'];
% bio_file=['roms_river_auto_new8_',num2str(year),'_yangtze_fennel_po4_8rivers.nc'];

tsuv_file=['roms_river_auto_new8_',num2str(year),'_yangtze80s_po4_8rivers.nc'];
bio_file=['roms_river_auto_new8_',num2str(year),'_yangtze80s_fennel_po4_8rivers.nc'];

% tsuv_file=['roms_river_auto_new8_',num2str(year),'_inYS_po4_8rivers.nc'];
% bio_file=['roms_river_auto_new8_',num2str(year),'_inYS_fennel_po4_8rivers.nc'];


nr=netcdf(bio_file,'r');
nc=netcdf(tsuv_file,'w');

temp=nc{'river_temp'}(:);
[time,s_rho,river]=size(temp);

theVarname = 'river_NH4';
nc{theVarname} = ncfloat('time','s_rho','river');
nc{theVarname}.long_name = ncchar('river runoff NH4');
nc{theVarname}.units = ncchar('millie Mole N meter-3');
nc{theVarname}.field = ncchar('river_NH4, scalar, series');
nc{theVarname}.time = ncchar('river_time');
nc{theVarname}(:)=nr{theVarname}(:);

theVarname = 'river_NO3';
nc{theVarname} = ncfloat('time','s_rho','river');
nc{theVarname}.long_name = ncchar('river runoff NO3');
nc{theVarname}.units = ncchar('millie Mole N meter-3');
nc{theVarname}.field = ncchar('river_NO3, scalar, series');
nc{theVarname}.time = ncchar('river_time');
nc{theVarname}(:)=nr{theVarname}(:);

theVarname = 'river_detritus';
nc{theVarname} = ncfloat('time','s_rho','river');
nc{theVarname}.long_name = ncchar('river runoff detritus');
nc{theVarname}.units = ncchar('millie Mole N meter-3');
nc{theVarname}.field = ncchar('river_detritus, scalar, series');
nc{theVarname}.time = ncchar('river_time');
nc{theVarname}(:)=nr{theVarname}(:);

theVarname = 'river_Oxyg';
nc{theVarname} = ncfloat('time','s_rho','river');
nc{theVarname}.long_name = ncchar('river runoff oxygen');
nc{theVarname}.units = ncchar('microMole L-1');
nc{theVarname}.field = ncchar('river_oxygen, scalar, series');
nc{theVarname}.time = ncchar('river_time');
nc{theVarname}(:)=nr{theVarname}(:);

theVarname = 'river_tPO4';
nc{theVarname} = ncfloat('time','s_rho','river');
nc{theVarname}.long_name = ncchar('river runoff PO4');
nc{theVarname}.units = ncchar('microMole L-1');
nc{theVarname}.field = ncchar('river_PO4, scalar, series');
nc{theVarname}.time = ncchar('river_time');
nc{theVarname}(:)=nr{theVarname}(:);

close(nc);close(nr);