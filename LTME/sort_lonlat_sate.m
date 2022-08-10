
close all
clear all
clc

np=netcdf('f:\avhrr\2006\AVHRR\amsr-avhrr-v2.20060101.nc');
lonr=np{'lon'}(:);
latr=np{'lat'}(:);
close(np)
st_lor = find(lonr <= 123.2 & lonr >= 123);
ed_lor = find(lonr <= 133.2 & lonr >= 133);
st_lar = find(latr <=  30.2 & latr >=  30);
ed_lar = find(latr <=  40.2 & latr >=  40);
lonr=lonr(st_lor:ed_lor);
latr=latr(st_lar:ed_lar);
[lonr latr] = meshgrid(lonr,latr);

[x y] = size(lonr);

lonr2 = reshape(lonr,1,x*y);
latr2 = reshape(latr,1,x*y);
data = [lonr2; latr2];
fid = fopen('position_sate.dat','w+');
fprintf(fid,'%12.7f%12.7f\n',data);
fclose(fid);