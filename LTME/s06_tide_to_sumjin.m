clc; clear all; 
close all

%%
fname = '../../\01_2dm(basemap)\grid\Sumjin_roms_v11.2dm';
% fname = 'data\GY_v1.2dm';
[C a1 a2 a3 a4  ] = textread(fname,'%s  %f %f %f %f %*[^\n]', 'headerlines',2); % 보통 *.2dm 파일은 headerlines 이 2줄임

data = [a1 a2 a3 a4];

cc = char(C);
for i = 1:length(C)
    if(cc(i,:) == 'ND ')
        % if(cc(i,:) == 'ND')
        idcc(i,:) = 1;
    end
end
ia = find(idcc ==1);
sub_grid = data(min(ia):max(ia),:);
lm = a2(1) -2 ;
ln =  length(sub_grid) / lm ;

% lon = sub_grid(:,2);
% lat = sub_grid(:,3);

X = reshape(sub_grid(:,2),lm,ln)';
Y = reshape(sub_grid(:,3),lm,ln)';

%--- 큰 모델에서 가져올 영역 선택 ----
gridfile='data\roms_grd.nc';
nt=netcdf(gridfile,'r');
m_lon = 335;
tlon=nt{'lon_rho'}(:); tlon = tlon(292:313, 301:m_lon);
tlat=nt{'lat_rho'}(:); tlat = tlat(292:313, 301:m_lon);
x = tlon;
y = tlat;

gridfile='data\roms_tide.nc';
nt=netcdf(gridfile,'r');
tide_Cangle=nt{'tide_Cangle'}(:); tide_Cangle = tide_Cangle(:,292:313, 301:m_lon);
tide_Cmax = nt{'tide_Cmax'}(:); tide_Cmax = tide_Cmax(:,292:313, 301:m_lon);
tide_Cmin = nt{'tide_Cmin'}(:); tide_Cmin = tide_Cmin(:,292:313, 301:m_lon);
tide_Cphase = nt{'tide_Cphase'}(:); tide_Cphase = tide_Cphase(:,292:313, 301:m_lon);
tide_Eamp = nt{'tide_Eamp'}(:); tide_Eamp = tide_Eamp(:,292:313, 301:m_lon);
tide_Ephase = nt{'tide_Ephase'}(:); tide_Ephase = tide_Ephase(:,292:313, 301:m_lon);
tide_period = nt{'tide_period'}(:); 


for i = 1:8
    z = squeeze(tide_Cangle(i,:,:));
    %------------------------------------
    Z(i,:,:)=griddata(x,y,z,X,Y); 
end
tide_Cangle = Z;

for i = 1:8
    z = squeeze(tide_Cmax(i,:,:));
    z(isnan(z)==0)= 0;
    %------------------------------------
    Z(i,:,:)=griddata(x,y,z,X,Y); 
end
tide_Cmax = Z;

for i = 1:8
    z = squeeze(tide_Cmin(i,:,:));
    z(isnan(z)==0)= 0;
    %------------------------------------
    Z(i,:,:)=griddata(x,y,z,X,Y); 
end
tide_Cmin = Z;

for i = 1:8
    z = squeeze(tide_Cphase(i,:,:));
    %------------------------------------
    Z(i,:,:)=griddata(x,y,z,X,Y); 
end
tide_Cphase = Z;
pcolor(squeeze(Z(2,:,:)));shading flat;colorbar;
for i = 1:8
    z = squeeze(tide_Eamp(i,:,:));
    
    if i == 1
        z = z+z*0.35;
    elseif i == 2
        z = z+z*0.33; % 0.15
    elseif i == 3
        z = z+z*0.35;
    elseif i ==4
        z = z+z*0.15;        
   elseif i == 5
       z = z+z*0.16;
    elseif i ==6
        z = z+z*0.11;
%     elseif i == 7
%         z = z-z*0.02;        
    elseif i == 8
        z = z-z*0.22;        
    end
    %------------------------------------
    Z(i,:,:)=griddata(x,y,z,X,Y); 
end
tide_Eamp = Z;

for i = 1:8
    z = squeeze(tide_Ephase(i,:,:));
    %------------------------------------
    Z(i,:,:)=griddata(x,y,z,X,Y); 
end
tide_Ephase = Z;
pcolor(squeeze(Z(1,:,:)));shading flat;colorbar;
% pcolor(squeeze(Z(1,:,:)));colorbar;
%% nc 파일로 만들기
nc=netcdf('sumjin_v11_tide_larged3.nc','clobber');

%--- Global attributes:
nc.date = ncchar('August_07_2016');
nc.start_tide_mjd = ncdouble(48622);
nc.components = ncchar('M2 S2 N2 K2 K1 O1 P1 Q1');
nc.type = ncchar('ROMS - Tidal forcing (sumjin_v11)');
% nc.grd_file = ncchar('H:\TEST\yw3km_orth3\grd\roms_grd_yw36orth3c.nc');
nc.author = ncchar('silver in SNU');

%--- Dimensions:
ncdim('tide_period',8,nc);
ncdim('eta_rho',176,nc); % 모델grid에 맞게 수정하기
ncdim('xi_rho',252,nc);

%--- Variables and attributes:
nc{'tide_period'} = ncdouble('tide_period'); %% 8 elements.
nc{'tide_period'}.long_name = ncchar('Tide angular periode');
nc{'tide_period'}.units = ncchar('Hour');
nc{'tide_period'}(:)=tide_period;

nc{'tide_Ephase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho'); %% 1506944 elements.
nc{'tide_Ephase'}.long_name = ncchar('Tidal elevation phase angle');
nc{'tide_Ephase'}.units = ncchar('Degrees');
nc{'tide_Ephase'}(:)=tide_Ephase;

nc{'tide_Eamp'} = ncdouble('tide_period', 'eta_rho', 'xi_rho'); %% 1506944 elements.
nc{'tide_Eamp'}.long_name = ncchar('Tidal elevation amplitude');
nc{'tide_Eamp'}.units = ncchar('Meter');
nc{'tide_Eamp'}(:) = tide_Eamp;

nc{'tide_Cmin'} = ncdouble('tide_period', 'eta_rho', 'xi_rho'); %% 1506944 elements.
nc{'tide_Cmin'}.long_name = ncchar('Tidal current ellipse semi-minor axis');
nc{'tide_Cmin'}.units = ncchar('Meter second-1');
nc{'tide_Cmin'}(:) = tide_Cmin;

nc{'tide_Cmax'} = ncdouble('tide_period', 'eta_rho', 'xi_rho'); %% 1506944 elements.
nc{'tide_Cmax'}.long_name = ncchar('Tidal current, ellipse semi-major axis');
nc{'tide_Cmax'}.units = ncchar('Meter second-1');
nc{'tide_Cmax'}(:) = tide_Cmax;
 
nc{'tide_Cangle'} = ncdouble('tide_period', 'eta_rho', 'xi_rho'); %% 1506944 elements.
nc{'tide_Cangle'}.long_name = ncchar('Tidal current inclination angle');
nc{'tide_Cangle'}.units = ncchar('Degrees between semi-major axis and East');
nc{'tide_Cangle'}(:) = tide_Cangle;
 
nc{'tide_Cphase'} = ncdouble('tide_period', 'eta_rho', 'xi_rho'); %% 1506944 elements.
nc{'tide_Cphase'}.long_name = ncchar('Tidal current phase angle');
nc{'tide_Cphase'}.units = ncchar('Degrees');
nc{'tide_Cphase'}(:) = tide_Cphase;

endef(nc)
close(nc)

%}





