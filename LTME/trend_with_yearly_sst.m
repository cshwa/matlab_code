clear all;clc;close all;

%read reference setting file
filepath ='E:\11.사업\장기생태_3단계\Data\위성자료/yearly/';
name1 = 'avhrr_';
file1 = strcat(filepath, 'avhrr_1982.nc');
nc = netcdf(file1);
row_lat = nc{'lat'}(:); %unit: degrees_E  
row_lon = nc{'long'}(:);%unit: degrees_N 
row_temp= nc{'temp'}(:);%unit: Celsius
% clear nc;
st_y = 1984;
en_y = 2015;
k = 1;
%creating nc file : yearly temp
for i = st_y:en_y
    temp_sum = zeros(size(row_temp));
    name2 = num2str(i);
        file_name = strcat(filepath,name1, name2,'.nc');
        nc = netcdf(file_name);
        temp_sum = temp_sum + nc{'temp'}(:);
        data = nc{'temp'}(:);
        to_ye_temp(:,:,k) = data(488:502,501:520);
        k = k+1;
        clear nc;
end

save to_ye_temp to_ye_temp







