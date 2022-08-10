% Calculate total kinetic energy
% units [m2/s2]  == energy per unit mass [J/kg]
% edited by YK.Kang

clc; clear all; close all;
domainname = 'MICT';
expname = 'MICT_1_30';
testname = 'test04'; 

gridpath = ['D:\ROMS\data\MICT\grid\1_30\',testname,'\'];
gridname = ['grid_',expname,'_',testname,'.nc'];
gname = [gridpath,gridname];
year = 1993; yy = num2str(year)

out_path = ['H:\',expname,'\KE\',testname,'\'];
        
g = grd([expname,'_',testname]);
h = g.h;
hc = g.hc;
theta_s = g.theta_s;
theta_b = g.theta_b;
N = g.N;
mask_rho = g.mask_rho ./ g.mask_rho;
lon_rho = g.lon_rho;
lat_rho = g.lat_rho;
angle = g.angle;
pm = g.pm;
pn = g.pn;
domaxis = [min(min(lon_rho)) max(max(lon_rho)) min(min(lat_rho)) max(max(lat_rho))];
[L, M] = size(lon_rho);
h = h .* mask_rho;

% index = find(lat_rho >= 40) ;
% load_path = 'D:\ROMS\eastsea\data\output\new_test9\2012\KE\';
% load_name = 'new_test9_TKE.mat';
% load([load_path, load_name])
index = find(lon_rho>=143);

for spin = 2
    spinup_name = [num2str(spin),'yr'];
    
    for month = 1 : 12
        mm = num2char(month,2);
        disp([spinup_name,'.',mm])
        filepath = ['H:\',expname,'\output\',testname,'\spinup\EnOI\',spinup_name,'\'];
        filename = ['\monthly_',yy,mm,'.nc'];
        
        fname = [filepath, filename];
        nc = netcdf(fname);
        u = nc{'u'}(:); v = nc{'v'}(:);
        u = u2rho_3d(u); v = v2rho_3d(v);
        u(u>100) = NaN; v(v>100) = NaN;
        u(index) = NaN; v(index) = NaN;
%         for l = 1 : N
%             uu = squeeze(u(l,:,:));
%             uu(index) = NaN;
%             u(l,:,:) = uu;
%             vv = squeeze(v(l,:,:));
%             vv(index) = NaN;
%             v(l,:,:) = vv;
%         end
        
%         h(index) = NaN;
        
        zeta = nc{'zeta'}(:);
        Hz = diff(zlevs(h,zeta,theta_s,theta_b,g.hc,N,'w'));
        
        % Kinetic energy = 1/2 * m * v2
        dx = 1./pm;
        dy = 1./pn;
        volume = dx .* dy .* (h+zeta) .* (mask_rho ./ mask_rho); % 면적 * 각 층의 수심(두께)
        total_volume = nansum(nansum(volume));
        
        for   k = 1 : N
            for j = 1 : L
                for i = 1 : M
                    if j == 1 | i == 1 | j == L | i == M
                        v2(k,j,i) = NaN;
                        
                    else
                        v2(k,j,i) = ((u(k,j,i) .* u(k,j,i)) + ...
                            (u(k,j,i+1) .* u(k,j,i+1)) +...
                            (v(k,j,i) .* v(k,j,i)) + ...
                            (v(k,j+1,i) .* v(k,j+1,i)));
                    end
                end
            end
        end
        
        for k = 1 : N
            m(k,:,:) = dx .* dy .* squeeze(Hz(k,:,:));
        end
        
        kinetic = 1/4 .* m .* v2;
        TKE(month,spin) = nansum(nansum(nansum(kinetic))) / total_volume;
    end
end
save([out_path, testname,'_TKE.mat'],'kinetic','TKE')


