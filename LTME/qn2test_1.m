% Creat NetCDFfile by using wetcdf toolbox
% Input file = 'quickscat.mat'
% Programmed by JongJin Park

load quikscat.mat

[wcid_n,status]=wccreate('test.nc','clobber');
xt = lon;
yt = lat;
tt = time;

len_x = length(xt);
len_y = length(yt);
len_t = length(tt);
[dimid_x,status]=wcdimdef(wcid_n,'X',len_x);
[dimid_y,status]=wcdimdef(wcid_n,'Y',len_y);
[dimid_t,status]=wcdimdef(wcid_n,'T',len_t);

[varid_xn,status]=wcvardef(wcid_n,'X','float',1,dimid_x);
[status,len_att] = size('longitude');
status = wcattput(wcid_n,varid_xn,'long_name','char',len_att,'longitude');
[status,len_att] = size('degrees');
status = wcattput(wcid_n,varid_xn,'units','char',len_att,'degrees');

[varid_yn,status]=wcvardef(wcid_n,'Y','float',1,dimid_y);
[status,len_att] = size('latitude');
status = wcattput(wcid_n,varid_yn,'long_name','char',len_att,'latitude');
[status,len_att] = size('degrees');
status = wcattput(wcid_n,varid_yn,'units','char',len_att,'degrees');

[varid_tn,status]=wcvardef(wcid_n,'T','float',1,dimid_t);
[status,len_att] = size('Time');
status = wcattput(wcid_n,varid_tn,'long_name','char',len_att,'Time');
[status,len_att] = size('days');
status = wcattput(wcid_n,varid_tn,'units','char',len_att,'days');

status=wcvarput(wcid_n,varid_xn,[0],[len_x],xt);
status=wcvarput(wcid_n,varid_yn,[0],[len_y],yt);
status=wcvarput(wcid_n,varid_tn,[0],[len_t],tt);

% eta
[varid_u,status]=wcvardef(wcid_n,'taux','float',3,[dimid_t,dimid_y,dimid_x]);
attvalue='N/m2';
[status,len_att] = size(attvalue);
status = wcattput(wcid_n,varid_u,'units','char',len_att,attvalue);
status = wcattput(wcid_n,varid_u,'missing_value','float',1,-1e34);
status = wcvarput(wcid_n,varid_u,[0,0,0],[len_t,len_y,len_x],u);

[varid_v,status]=wcvardef(wcid_n,'tauy','float',3,[dimid_t,dimid_y,dimid_x]);
attvalue='N/m2';
[status,len_att] = size(attvalue);
status = wcattput(wcid_n,varid_v,'units','char',len_att,attvalue);
status = wcattput(wcid_n,varid_v,'missing_value','float',1,-1e34);
status = wcvarput(wcid_n,varid_v,[0,0,0],[len_t,len_y,len_x],v);

wcclose(wcid_n);
disp('done')
