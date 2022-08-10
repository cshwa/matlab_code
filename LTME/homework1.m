clear all
%u=ncread('analysis.mon.mean.nc','u');
v=ncread('vwnd.mon.mean.nc','vwnd');
w=ncread('omega.mon.mean.nc','omega');
%zonalw=ncread('analysis.mon.mean.nc','Zw');
%zv=ncread('analysis.mon.mean.nc','Zv');
%llat=ncread('omega.mon.mean.nc','lat');
%zeta=ncread('analysis.mon.mean.nc','zeta');
p=ncread('analysis.mon.mean.nc','level');
%omegap=ncread('omega.mon.mean.nc','level');
%phi=ncread('hgt.mon.mean.nc','hgt');
% pressure units -> milibar.
% omega units -> pascal/sec.  101325 pascal=1atm.
attvalue = ncreadatt('omega.mon.mean.nc','omega','units')
attvalue = ncreadatt('omega.mon.mean.nc','level','units')
attvalue = ncreadatt('vwnd.mon.mean.nc','level','units')

pascal=101325;
p=p * (1/1013.25) * (101325/1);
% for i=1:17
%     dp(i)=p(18-i)
% end
% for i=1:17
%     p(i)=dp(i)
% end
zv(1:73,1:17,1:432)=0;
zw(1:73,1:17,1:432)=0;
zphi(1:73,1:17,1:432)=0;
for j=1:73
    for k=1:17
        for l=1:432
            for i=1:144
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


%for i=1:432
%   tzw(:,:,i)=zw(:,:,i).';
%    tzv(:,:,i)=zv(:,:,i).';
%    tzeta(:,:,i)=zeta(:,:,i).';
%end
%p=p.'
% pi=3.14159265358979323846
% omega=0.000072921
% for j=1:73
%     lat(j)=llat(74-j);
% end
% for j=1:73
%     f(j)=2.0*omega*sin(lat(j)*pi/180.0);
% end


dy=277987.3
%dy=25000
for j=1:72
    for k=1:16
        for l=6:8
            zeta(j,k,l)=(zw(j+1,k,l)-zw(j,k,l))/(dy) - (zv(j,k+1,l)-zv(j,k,l))/(p(k+1)-p(k));
        end
    end
end
for j=2:72
    for k=2:16
        for l=6:8
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
for n=1:73
    for i=0:16;
        tp(73*i+n)=p(i+1);
    end
end
for i=1:1241
    for n=1:1241
        ttp(n,i)=tp(i);
    end
end

for i=74:1168
   chi(i,i)= -(2/((dy)^2) + 2/((ttp(i,i+73)-ttp(i,i-73))*(ttp(i,i+73)-ttp(i,i))) + 2/((ttp(i,i+73)-ttp(i,i-73))*(ttp(i,i)-ttp(i,i-73))));
   chi(i+1,i) = 2/((ttp(i,i+73)-ttp(i,i-73))*(ttp(i,i)-ttp(i,i-73)));
   chi(i,i+1) = 2/((ttp(i,i+73)-ttp(i,i-73))*(ttp(i,i+73)-ttp(i,i)));
   chi(73+i,i)= 1/(dy^2) ;
   chi(i,17+i)= 1/(dy^2);
end
for i=1:73
   chi(i,i)= -(2/((dy)^2) + 1/((ttp(i,i+73)-ttp(i,i))*(ttp(i,i+73)-ttp(i,i))));
   chi(i+1,i) = 0;
   chi(i,i+1) = 1/((ttp(i,i+73)-ttp(i,i))*(ttp(i,i+73)-ttp(i,i)));
   chi(73+i,i)= 1/(dy^2) ;
   chi(i,17+i)= 1/(dy^2);
end
for i=1169:1241
   chi(i,i)= -(2/((dy)^2) + 1/((ttp(i,i)-ttp(i,i-73))*(ttp(i,i)-ttp(i,i-73))));
    if i<1241
      chi(i+1,i) = 1/((ttp(i,i)-ttp(i,i-73))*(ttp(i,i)-ttp(i,i-73)));
      chi(i,i+1) = 0;
    end
    if i<1225
        chi(i,17+i)= 1/(dy^2);
    end 
end

%%%%%%%%%%%%%%%%%%%%%%%%
%zeta=zeta.*100000;
%%%%%%%%%%%%%%%%%%%%%%%%%

% dp=1000
% for i=74:1168
%    chi(i,i)= -2/dy-2/dp;
%    chi(i+1,i) = 1/dp;
%    chi(i,i+1) = 1/dp;
%    chi(73+i,i)=1/dy ;
%    chi(i,17+i)= 1/dy;
% end
% for i=1:73
%    chi(i,i)= -2/dy-1/dp;
%    chi(i+1,i) = 0;
%    chi(i,i+1) = 1/dp;
%    chi(73+i,i)= 1/dy ;
%    chi(i,17+i)= 1/dy;
% end
% for i=1169:1241
%    chi(i,i)= -2/dy-1/dp;
%     if i<1241
%       chi(i+1,i) = 1/dp;
%       chi(i,i+1) = 0;
%     end
%     if i<1225
%         chi(i,17+i)= 1/dy;
%     end 
% end




for n=1:16
    chi(73*n+1,17*n) =0;
    chi(73*n,17*n+1) =0;
end
invchi=inv(chi);
for t=6:8
    for j=0:16
        for i=1:73
            tzeta(j*73+i,t)=zeta(i,j+1,t);
        end
    end
end
% for t=6:8
%     for i=1:1241
%         ttzeta(i,i,t)=tzeta(i,t); 
%     end
% end
for t=6:8
    chichi(:,t)=invchi*tzeta(:,t);
end
for t=6:8
    for i=1:73
        for j=0:16
            finalchi(i,j+1,t)=chichi(j*73+1,t);
           % if finalchi(i,j+1,t)>150
          %      finalchi(i,j+1,t)=0;
        %    elseif finalchi(i,j+1,t)<-150
         %       finalchi(i,j+1,t)=0;
       %     end
        end
    end
end
% for j=1:73
%     for k=1:17
%         phif(j,k)=zphi(j,k,50)/f(j);
%     end
% end
%surf(phif);
finalchi=finalchi/8
surf(zeta(:,2:11,8))
figure;
surf(finalchi(:,2:11,8))


