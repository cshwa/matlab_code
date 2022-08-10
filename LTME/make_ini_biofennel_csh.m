clc;clear all;close all

inifile='pohang_ini_30_fix_1m_1m_1c_v11_po4_fill_v2.nc';

nw=netcdf(inifile,'write');

temp=nw{'temp'}(:);
[z,y,x]=size(temp);

nw('s_rho')=z;
nw('eta_rho')=y;
nw('xi_rho')=x;
nw('ocean_time')=1;
nw{'ocean_time'}(:)=864000;
% nw{'ocean_time'}(:)=12873600; %% 149 = 5/29
% nw{'ocean_time'}(:)=19526400; %% 226 = 8/14

%%%%%  - initial 에 2-layer를 주기 위해 NO3,NH4는 make_ini_pohang.m 에서 넣어줌. %%%%%
% nw{'NO3'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nw{'NO3'}.long_name = ncchar('Nitrate');
% nw{'NO3'}.units = ncchar('mMol N m-3');
% nw{'NO3'}.fields = ncchar('NO3, scalar, series');
% nw{'NO3'}(:)=13.5;

% nw{'NH4'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nw{'NH4'}.long_name = ncchar('Ammonium');
% nw{'NH4'}.units = ncchar('mMol N m-3');
% nw{'NH4'}.fields = ncchar('NH4, scalar, series');
% nw{'NH4'}(:)=160;

% nw{'tPO4'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nw{'tPO4'}.long_name = ncchar('Phosphate');
% nw{'tPO4'}.units = ncchar('milimole PO4 m-3');
% nw{'tPO4'}.fields = ncchar('tPO4, scalar, series');
% nw{'tPO4'}(:)=0.2;
% 
% nw{'chlorophyll'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nw{'chlorophyll'}.long_name = ncchar('chlorophyll');
% nw{'chlorophyll'}.units = ncchar('mg C m-3');
% nw{'chlorophyll'}.fields = ncchar('CHLA, scalar, series');
% nw{'chlorophyll'}(:)=0.5;

nw{'phytoplankton'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho');
nw{'phytoplankton'}.long_name = ncchar('phytoplankton');
nw{'phytoplankton'}.units = ncchar('mMol N m-3');
nw{'phytoplankton'}.fields = ncchar('PHYTO, scalar, series');
nw{'phytoplankton'}(:)=0.05;

nw{'zooplankton'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nw{'zooplankton'}.long_name = ncchar('zooplankton');
nw{'zooplankton'}.units = ncchar('mMol N m-3');
nw{'zooplankton'}.fields = ncchar('ZOO, scalar, series');
nw{'zooplankton'}(:)=0.01;

nw{'LdetritusN'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nw{'LdetritusN'}.long_name = ncchar('large fraction nitrogen detritus concentration');
nw{'LdetritusN'}.units = ncchar('mMol N m-3');
nw{'LdetritusN'}.fields = ncchar('LdetritusN, scalar, series');
nw{'LdetritusN'}(:)=0;

nw{'SdetritusN'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nw{'SdetritusN'}.long_name = ncchar('small fraction nitrogen detritus concentration');
nw{'SdetritusN'}.units = ncchar('mMol N m-3');
nw{'SdetritusN'}.fields = ncchar('SdetritusN, scalar, series');
nw{'SdetritusN'}(:)=0;

nw{'LdetritusC'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nw{'LdetritusC'}.long_name = ncchar('large fraction carbon detritus concentration');
nw{'LdetritusC'}.units = ncchar('mMol C m-3');
nw{'LdetritusC'}.fields = ncchar('LdetritusC, scalar, series');
nw{'LdetritusC'}(:)=0;

nw{'SdetritusC'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nw{'SdetritusC'}.long_name = ncchar('small fraction small detritus concentration');
nw{'SdetritusC'}.units = ncchar('mMol C m-3');
nw{'SdetritusC'}.fields = ncchar('SdetritusC, scalar, series');
nw{'SdetritusC'}(:)=0;

nw{'TIC'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nw{'TIC'}.long_name = ncchar('total inorganic carbon');
nw{'TIC'}.units = ncchar('mMol Carobn m-3');
nw{'TIC'}.fields = ncchar('TIC, scalar, series');
nw{'TIC'}(:)=1;

nw{'alkalinity'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
nw{'alkalinity'}.long_name = ncchar('total alkalinity');
nw{'alkalinity'}.units = ncchar('milliequivalents meter-3');
nw{'alkalinity'}.fields = ncchar('alkalinity, scalar, series');
nw{'alkalinity'}(:)=2;

% nw{'oxygen'} = ncdouble('ocean_time', 's_rho', 'eta_rho', 'xi_rho'); 
% nw{'oxygen'}.long_name = ncchar('dissolved oxygen concentration');
% nw{'oxygen'}.units = ncchar('millimole_oxygen meter-3');
% nw{'oxygen'}.fields = ncchar('oxygen, scalar, series');
% nw{'oxygen'}(:)=10;

close(nw);
