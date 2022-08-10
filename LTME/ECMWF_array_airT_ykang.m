% clc; clear all; close all;
% var = 'v10';
% varname = 'v10'
% 
% for y = 2017
%     num = 1;
%     year = num2str(y); % year
%     fpath = ['D:\ROMS\eastsea\data\air_forcing\ECMWF\',var,'\'];
%     fname = [fpath, 'ECMWF_Interim_',var,'_',year,'_re.nc'];
%     create_ECMWF(fname,y,var,varname)
%     nk = netcdf(fname, 'w');
%     
%     %     for m = 1:12;
%     %         mon = num2char(m,2); %month
%     %         disp(['year is ', year]);
%     %         disp(['month is ', mon]);
%     
%     filename = [fpath,'ECMWF_Interim_',var,'_',year,'.nc'];
%     nc = netcdf(filename);
%     time = nc{'time'}(:);
%     longitude = nc{'longitude'}(:); latitude = nc{'latitude'}(:);
%     temp_scale_factor = nc{varname}.scale_factor(:);
%     temp_add_offset = nc{varname}.add_offset(:);
%     date_num = datenum(1900,1,1) + time/24;
%     date = datevec(date_num);
%     yy = date(:,1); mm = date(:,2);
%     day = date(:,3); hour = date(:,4);
%     
%     for m = 1:12
%         mon = num2char(m,2)
%         
%         for dd = 1:eomday(y,m)
%             for h = 0:6:18
%                 idx_day = find(yy == y & mm == m & day == dd & hour == h);
%                 if length(idx_day) ~= 0
%                     t2m_o = nc{varname}(idx_day,:,:);
%                     nk{'time'}(num,1) = time(idx_day);
%                     check_time(num,1) = time(idx_day);
%                     nk{varname}(num,:,:) = t2m_o.*temp_scale_factor+temp_add_offset;
%                     num = num+1;
%                 end
%             end
%         end
%         
%         clear vars idx_day  m h
%     end
%     
%     for m = 1:12
%         mon = num2char(m,2) %month
%         %         disp(['year is ', year]);
%         %         disp(['month is ', mon]);
%         %         filepath = ['D:\ROMS\eastsea\forcing\air\','ECMWF_Interim_airT_',year,'.nc'];
%         %         nc = netcdf(filepath);
%         %         time = nc{'time'}(:);
%         %         longitude = nc{'longitude'}(:); latitude = nc{'latitude'}(:);
%         %         temp_scale_factor = nc{'t2m'}.scale_factor(:);
%         %         temp_add_offset = nc{'t2m'}.add_offset(:);
%         %         date_num = datenum(1900,1,1) + time/24;
%         %         date = datevec(date_num);
%         %         yy = date(:,1); mm = date(:,2);
%         %         day = date(:,3); hour = date(:,4);
%         %
%         for dd = 1:eomday(y,m)
%             for h = 3:6:21
%                 idx_day = find(yy == y & mm == m & day == dd & hour == h);
%                 if length(idx_day) ~= 0
%                     t2m_o = nc{varname}(idx_day,:,:);
%                     nk{'time'}(num,1) = time(idx_day);
%                     check_time(num,1) = time(idx_day);
%                     nk{varname}(num,:,:) = t2m_o.*temp_scale_factor+temp_add_offset;
%                     num = num+1;
%                 end
%             end
%         end
%     end
%     
%     %
%     disp(['end']);
%     nk{varname}(num,:,:) = nc{varname}(end,:,:).*temp_scale_factor+temp_add_offset;
%     nk{'time'}(num,1) = time(end);
%     check_time(num,1) = time(end);
%     nk{'longitude'}(:) = longitude;
%     nk{'latitude'}(:) = latitude;
%     
%     close(nc)
%     close(nk)
%     clc; clear all;
%     
% end

%% make ssrd file
clc; clear all; close all;
var = 'ssrd';
varname = 'ssrd'
grdfile='D:\ROMS\eastsea\data\grid\2012\new_test9\new_etopo1_Eastsea_test9.nc';
a = 1;

for y = 2009
    disp(['year = ',num2str(y)])
    num = 1;
    year = num2str(y); % year
    
    ly = (( mod(y,4) == 0 & mod(y,100) ~= 0 ) | mod(y,400) == 0);
    if ly,    days = 366;
    else     days = 365;
    end
    
    fpath = ['D:\ROMS\eastsea\data\air_forcing\ECMWF\',var,'\'];
    %     fname = [fpath, 'ECMWF_Interim_',var,'_',year,'.nc'];
    %     create_ECMWF(fname,y,var,varname)
    %     nk = netcdf(fname, 'w');
    
    filename = [fpath,'ECMWF_Interim_',var,'_',year,'.nc'];
    nc = netcdf(filename);
    time = nc{'time'}(:);
    longitude = nc{'longitude'}(:); latitude = nc{'latitude'}(:);
    temp_scale_factor = nc{varname}.scale_factor(:);
    temp_add_offset = nc{varname}.add_offset(:);
    date_num = datenum(1900,1,1) + time/24;
    date = datevec(date_num);
    yy = date(:,1); mm = date(:,2);
    day = date(:,3); hour = date(:,4);
    
    for m = 1:12
        mon = num2char(m,2)
        
        for dd = 1:eomday(y,m)
            index_12 = find(yy == y & mm == m & day == dd & hour == 12);
            
            if m == 12 && dd == eomday(y,m)
                index_24 = find(yy == y+1 & mm == 1 & day == 1 & hour == 0);
            elseif dd == eomday(y,m)
                %                 index_24 = find(yy == y & mm == m+1 & day == 1 & hour == 0);
                index_24 = length(time);
            else
                index_24 = find(yy == y & mm == m & day == dd+1 & hour == 0);
            end
            
            var_12 = nc{varname}(index_12,:,:);
            var_12 = var_12 .* temp_scale_factor+temp_add_offset;
            var_24 = nc{varname}(index_24,:,:);
            var_24 = var_24 .* temp_scale_factor+temp_add_offset;
            daily(a,:,:) = (var_12 + var_24) ./ 86400;
            a = a+1;
            
        end
    end


grd=netcdf(grdfile);
lat_rho=grd{'lat_rho'}(:);
lon_rho=grd{'lon_rho'}(:);
close(grd);
[X,Y]=meshgrid(longitude,latitude);

for k= 1 : 1 : days
    disp(['day = ',num2str(k)])
    Z=squeeze(daily(k,:,:));
    Zi=griddata(X, Y, Z, lon_rho, lat_rho);
    swrad_daily(k,:,:)=Zi;
end

save(['swrad_',year,'_daily.mat'],'swrad_daily')
end
close all; clear all; clc;
days=365;
time=[0.5:1:days-0.5];
load swrad_2009_daily.mat
fname=['ECMWF_post_2009_swrad.nc'];
create_roms_forcing_SWD(fname,swrad_daily,time,days,days)







