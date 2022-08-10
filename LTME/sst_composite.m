function [error]=sst_composite

jnum = 60;
inum = 60;
tstart = 3;
tnum = 1;
twin=2;

xdlen=10;
ydlen=10;
err_t=0.5;

time_units = ncreadatt('sample.nc','TIME','units');
time_origins = ncreadatt('sample.nc','TIME','time_origin');

nc = netcdf.create('test.nc','clobber');
dimid_time = netcdf.defDim(nc,'time',netcdf.getConstant('NC_UNLIMITED'));
dimid_lon = netcdf.defDim(nc,'lon',inum);
dimid_lat = netcdf.defDim(nc,'lat',jnum);

varid_lon = netcdf.defVar(nc,'lon','double',[dimid_lon]);
netcdf.putAtt(nc,varid_lon,'long_name','Longitude');
netcdf.putAtt(nc,varid_lon,'units','degree_e');
netcdf.putAtt(nc,varid_lon,'cartesian_axis','X');

varid_lat = netcdf.defVar(nc,'lat','double',[dimid_lat]);
netcdf.putAtt(nc,varid_lat,'long_name','Latitude');
netcdf.putAtt(nc,varid_lat,'units','degree_n');
netcdf.putAtt(nc,varid_lat,'cartesian_axis','Y');

varid_time = netcdf.defVar(nc,'time','double',[dimid_time]);
netcdf.putAtt(nc,varid_time,'long_name','Time');
netcdf.putAtt(nc,varid_time,'units',time_units);
netcdf.putAtt(nc,varid_time,'origins',time_origins);
netcdf.putAtt(nc,varid_time,'cartesian_axis','T');

varid_sst = netcdf.defVar(nc,'SST','float',[dimid_lon,dimid_lat,dimid_time]);
netcdf.putAtt(nc,varid_sst,'long_name','Observed SST');
netcdf.putAtt(nc,varid_sst,'units','degree C');
netcdf.putAtt(nc,varid_sst,'missing_value',-1.e34);

varid_ana = netcdf.defVar(nc,'Tana','float',[dimid_lon,dimid_lat,dimid_time]);
netcdf.putAtt(nc,varid_ana,'long_name','Analyzed SST');
netcdf.putAtt(nc,varid_ana,'units','degree C');
netcdf.putAtt(nc,varid_ana,'missing_value',-1.e34);

varid_err = netcdf.defVar(nc,'Terr','float',[dimid_lon,dimid_lat,dimid_time]);
netcdf.putAtt(nc,varid_err,'long_name','Analysys Error');
netcdf.putAtt(nc,varid_err,'units','degree C');
netcdf.putAtt(nc,varid_err,'missing_value',-1.e34);

lat = ncread('sample.nc','LAT',601-ydlen,jnum+2*ydlen);
lon = ncread('sample.nc','LON',1,inum+2*xdlen);
time = ncread('sample.nc','TIME',1,5);

temp = ncread('sample.nc','TEMPERATURE',[1, 601-ydlen, 1],[inum+2*xdlen,jnum+2*ydlen, 5]);

index = find(isnan(temp));

size(index)

temp(index) = -1.e34;

tana = ones(jnum,inum);
terr = ones(jnum,inum);

for i=xdlen+1:xdlen:inum+xdlen
    i=i
  for j=ydlen+1:ydlen:jnum+ydlen
    for m=tstart:tstart+tnum-1
        
        iini = max(i-xdlen,1);
        iend = i+2*xdlen-1;
        jini = max(j-10,1);
        jend = j+2*ydlen-1;
        tini = m-twin;
        tend = m+twin;
        ic = iend-iini+1;
        jc = jend-jini+1;
        mc = tend-tini+1;

        mlat = mean(lat(j:j+xdlen-1));
        mlon = mean(lon(i:i+ydlen-1));

        if(iini == 1 & ic == xdlen+10)
            idata = 1;
        else
            idata = 11;
        end
        
        x=dist(mlat(ones(1,ic)),lon(iini:iend))/1000.;
        xi = ones(1,ic);
        xi(1) = 0.0;
        for ii=2:ic
            xi(ii) = xi(ii-1)+x(ii-1);
        end
        
        xi = xi-0.5*(xi(idata+xdlen/2-1)+xi(idata+xdlen/2));

        if(jini == 1 & jc == ydlen+10)
            jdata = 1;
        else
            jdata = 11;
        end

        y=dist(lat(jini:jend),mlon(ones(1,jc)))/1000.;
        yi = ones(jc,1);
        yi(1) = 0.0;
        for ii=2:jc
            yi(ii) = yi(ii-1)+y(ii-1);
        end
       
        yi = yi-0.5*(yi(jdata+ydlen/2-1)+yi(jdata+ydlen/2));

        ti = time(tini:tend);
        ti = ti-ti(twin+1);

        xdi = ones(ic,jc,mc);
        ydi = ones(ic,jc,mc);
        tdi = ones(ic,jc,mc);
        
        for it=1:ic;for jt=1:jc;for mt=1:mc
            xdi(it,jt,mt)=xi(it);
            ydi(it,jt,mt)=yi(jt);
            tdi(it,jt,mt)=ti(mt);
        end;end;end

        tpre = temp(iini:iend,jini:jend,tini:tend);
        
        index = find(tpre>-100.);
        x=xdi(index);
        y=ydi(index);
        t=tdi(index);
        z=tpre(index);

        for it=idata:idata+xdlen-1;for jt=jdata:jdata+ydlen-1
            xmi(it-idata+1,jt-jdata+1)=xi(it);
            ymi(it-idata+1,jt-jdata+1)=yi(jt);
        end;end

        [zi,err]=oi3d(x,y,t,z,xmi,ymi,50,50,5,err_t);
        tana(i-xdlen:i-1,j-ydlen:j-1)=zi;
        terr(i-xdlen:i-1,j-ydlen:j-1)=err;
    end
  end
end

index = find(abs(terr) > 0.05*err_t);
tana(index) = -1.e34;

netcdf.endDef(nc);

netcdf.putVar(nc,varid_lon,lon(xdlen+1:inum+xdlen));
netcdf.putVar(nc,varid_lat,lat(ydlen+1:jnum+ydlen));
netcdf.putVar(nc,varid_time,0,5,time);
netcdf.putVar(nc,varid_sst,temp(xdlen+1:inum+xdlen,ydlen+1:jnum+ydlen,1:5));
netcdf.putVar(nc,varid_ana,[0 0 2],[inum, jnum, 1], tana);
netcdf.putVar(nc,varid_err,[0 0 2],[inum, jnum, 1], terr);
netcdf.close(nc);

return;
