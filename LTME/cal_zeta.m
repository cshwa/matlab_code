clear all
 v=ncread('vwnd.mon.mean.nc','vwnd');
 w=ncread('omega.mon.mean.nc','omega');
p=ncread('omega.mon.mean.nc','level');
pascal=101325;
 p=p * (1/1013.25) * (101325/1);
%  for i=1:17
%     tp(i)=p(18-i);
%  end
%  for i=1:17
%      p(i)=tp(i);
%  end
zv(1:73,1:17,1:432)=0;
zw(1:73,1:17,1:432)=0;
%zphi(1:73,1:17,1:432)=0;
for j=1:73;
    for k=1:17;
        for l=1:432;
            for i=1:144;
                zv(j,k,l)=zv(j,k,l)+v(i,j,k,l);
                zw(j,k,l)=zw(j,k,l)+w(i,j,k,l);
%                zphi(j,k,l)=zphi(j,k,l)+phi(i,j,k,l);
            end
            zv(j,k,l)=zv(j,k,l)/144;
            zw(j,k,l)=zw(j,k,l)/144;
%            zphi(j,k,l)=zphi(j,k,l)/144;
        end
    end
end



dy=277987.3;

for j=2:72;
    for k=2:16;
        for l=1:432;
            zeta(j,k,l)=(zw(j+1,k,l)-zw(j-1,k,l))/(2*dy) - (zv(j,k+1,l)-zv(j,k-1,l))/(p(k+1)-p(k-1));
            %zeta(j,k,l)=-(zv(j,k+1,l)-zv(j,k,l))/(p(k+1)-p(k));
        end
    end
end
for j=2:72;
    for k=2:16;
        for l=1:432;
            if j==2
                zeta(j-1,k,l)=zeta(j,k,l);
            elseif j==72
                zeta(j+1,k,l)=zeta(j,k,l);
            end
            if k==2
                zeta(j,k-1,l)=zeta(j,k,l);
            elseif k==16
                zeta(j,k+1,l)=zeta(j,k,l);
            end
        end
    end
end

zeta8=zeta(:,:,8);



ncid= netcdf.create('reanalysis.mon.mean.nc','64BIT_OFFSET');
dimid=netcdf.defDim(ncid,'lat',73);
dimid2=netcdf.defDim(ncid,'level',17);
dimid3=netcdf.defDim(ncid,'time',0);

varid = netcdf.defVar(ncid,'zv','double',[dimid,dimid2,dimid3]);
varid2 = netcdf.defVar(ncid,'zw','double',[dimid,dimid2,dimid3]);
varid3 = netcdf.defVar(ncid,'zeta','double',[dimid,dimid2,dimid3]);

netcdf.endDef(ncid);

netcdf.putVar(ncid,varid,[0,0,0],[73,17,432],zv);
netcdf.putVar(ncid,varid2,[0,0,0],[73,17,432],zw);
netcdf.putVar(ncid,varid3,[0,0,0],[73,17,432],zeta);

netcdf.close(ncid);

ncdisp('reanalysis.mon.mean.nc');

zv=ncread('reanalysis.mon.mean.nc','zv');
zw=ncread('reanalysis.mon.mean.nc','zw');
zeta=ncread('reanalysis.mon.mean.nc','zeta');
dy=277987.3;
chi(1:73,1:17,1:432)=0;
C(1:73,1:17,1:432)=0;
D(1:73,1:17,1:432)=0;
E(1:73,1:17,1:432)=0;
chi2(1:73,1:17,1:432)=0;
dp=2000
for m=1:2
    for l=1:10
        for j=2:72
            for k=2:16
                   C(j,k,l) = (chi(j+1,k,l)-chi(j-1,k,l)) / (dy^2);
                   D(j,k,l) = chi(j,k+1,l) / ((p(k+1)-p(k-1))/2.) /(p(k+1)-p(k)) + chi(j,k-1,l) / ((p(k+1)-p(k-1))/2.) /(p(k)-p(k-1));
                   E(j,k,l) = (-2.)/(dy^2) + (1.)/((p(k+1)-p(k-1))/2.) /(p(k+1)-p(k)) + (1.)/((p(k+1)-p(k-1))/2.) /(p(k)-p(k-1));
%                  C(j,k,l) = (chi(j+1,k,l)-chi(j-1,k,l)) / (dy^2);
%                  D(j,k,l) = chi(j,k+1,l) / (dp) /(dp) + chi(j,k-1,l) / (dp) /(dp);
%                   E(j,k,l) = (-2.)/(dy^2) + (1.)/(dp) /(dp) + (1.)/(dp) /(dp);
                chi2(j,k,l) = (zeta(j,k,l)-C(j,k,l)-D(j,k,l))/E(j,k,l);
                if(j==2)
                    chi2(j-1,k,l)=chi2(j,k,l);
                end
                if(j==72)
                    chi2(j+1,k,l)=chi2(j,k,l);
                end
                if(k==2) 
                    chi2(j,k-1,l)=chi2(j,k,l);
                end
                if(k==16) 
                    chi2(j,k+1,l)=chi2(j,k,l);
                end
            end
        end
        for j=1:73
            for k=1:17
                for l=1:10
                  chi(j,k,l)=chi2(j,k,l);
                end
            end
        end
    end
end

% for l=1:3;
%     for j=1:73;
%         for k=1:17;
%             if(j==1)
%                 chi(j,k,l)=chi(2,k,l);
%             end
%             if(j==73)
%                 chi(j,k,l)=chi(73,k,l);
%             end
%             if(k==1) 
%                 chi(j,k,l)=chi(j,2,l);
%             end
%             if(k==17) 
%                 chi(j,k,l)=chi(j,16,l);
%             end
%         end
%     end
% end
C1=C(:,:,1);
D1=D(:,:,1);
E1=E(:,:,1);
zeta1=zeta(:,:,1);
chi1=chi(:,:,1);
chi2=chi(:,:,2);
chi3=chi(:,:,3);
surf(zeta(:,2:11,8))
figure;
surf(chi(:,2:11,8))