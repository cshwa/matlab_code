clc;clear all; close all;

data=load('timeseris_24h_ctd.txt');

max_dep=max(data(:,2));

n=0;
for i=1:1:length(data)-1
    if data(i,1)~=data(i+1,1) || i==length(data)-1
        max2_dep=max(data(i-n:i,2));
        gap_dep=max_dep-max2_dep;
        data(i-n:i,2)=data(i-n:i,2)+gap_dep;
        n=0;
    else
        n=n+1;
    end
end

% dep_cutdata=[data(1,2)];
% for itic=1:1:length(data)-1
%     if data(itic,1)==data(itic+1,1)
%         continue
%     else
%         dep_cutdata=[dep_cutdata,data(itic+1,2)];
%     end
% end


data=data(:,2:4);

asc_data=sortrows(data);
asc_data(:,1)=floor(asc_data(:,1));

leng=length(asc_data);
f_data=[];
k=0;n=1;
for i=1:1:(leng-1)
if asc_data(i,1)==asc_data(i+1,1) && i ~= leng -1
    n=n+1;
else if asc_data(i,1)~=asc_data(i+1,1) && asc_data(i,1)<asc_data(i+1,1)
    k=k+1;
    f_data(k,1)=asc_data(i,1);
    f_data(k,2)=sum(asc_data(i-n+1:i,2))/n;
    f_data(k,3)=sum(asc_data(i-n+1:i,3))/n;
    n=1;
    else if asc_data(i,1)~=asc_data(i+1,1) && asc_data(i,1)>asc_data(i+1,1)
        k=k+1;
    f_data(k,1)=asc_data(i,1);
    f_data(k,2)=sum(asc_data(i-n+1:i,2))/n;
    f_data(k,3)=sum(asc_data(i-n+1:i,3))/n;
        else if i == leng -1
             k=k+1;
    f_data(k,1)=asc_data(i,1);
    f_data(k,2)=sum(asc_data(i-n+1:i,2))/n;
    f_data(k,3)=sum(asc_data(i-n+1:i,3))/n;           
        break
            end
    end
end
end
end

%calculating density-------------------------------------------------------

tem=f_data(:,2);sal=f_data(:,3);
f_data(:,4)=999.842594+6.793952.*10.^(-2).*tem-9.095290*10^(-3).*tem.^(2)...
+1.001685*10^(-4).*tem.^(3)-1.120083*10^(-6).*tem.^(4)+6.536332*10^(-9).*tem.^(5)...
+8.24493*10^(-1).*sal-4.0899*10^(-3).*tem.*sal+7.6438*10^(-5).*tem.^2.*sal...
-8.2467*10^(-7).*tem.^3.*sal+5.3875*10^(-9).*tem.^(4).*sal-5.72466*10^(-3).*sal.^(3/2)...
+1.0227*10^(-4).*tem.*sal.^(3/2)-1.6546*10^(-6).*tem.^(2).*sal.^(3/2)+4.8314*10^(-4).*sal.^2; 
%--------------------------------------------------------------------------

plot(f_data(3:end,4),-[0:1:7]);