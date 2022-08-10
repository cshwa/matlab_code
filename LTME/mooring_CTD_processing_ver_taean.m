%------------------- 실험실 CTD rawdata 자료 처리  ------------------------%
% Column of data: Date Time Press Temp Cond Sal
%-------------------------------------------------------------------------%
% 1. input the number of ctd files
% 2. input longitude, latitude
% 3. file opan and text scan
clc;clear all;close all;

filenum=input('how many casts ?');
% filenum=str2double(filenum);

kk=0;

for fi=1:1:filenum
    [f, p]=uigetfile('*.*','select cast data file');
    filedir=[p,f];
    disp(filedir);
    %lat=input('latitude? : ');
    %lon=input('longitude? : ');
    time = input('Time? : ');
    %lat=str2num(lat);lon=str2num(lon);

    fid=fopen(filedir);
    c_header=textscan(fid,'%s %s %s %s %s %s ',1);
    c_data=textscan(fid,'%s %s %f %f %f %f' );
    fclose(fid);
    % c_data{3}=floor(c_data{3});
    data=[c_data{3},c_data{4},c_data{5},c_data{6}];
    maxdepth=max(data(:,1));
    check_dep=find(maxdepth == data(:,1));
    data=data(1:check_dep,:);
    data=sortrows(data,1);
    leng=length(data);
    
    k=0;
    n=0;
    out_zero=1;
    i1=0;
    
    % find mean between 0 <= depth < 0.5
    while out_zero == 1
     i1=i1+1;
        % use data only under the depth 0.5m 
        if data(i1,1) <= 0.5
            n=n+1;
        % if as the depth over 0.5m, find mean value 0m~ under 0.5m   
        else
            k=k+1;
            kk=kk+1;
            f_data(kk,1)=time;
            f_data(kk,2)=time;
            f_data(kk,3)=-floor(data(i1,1));
            f_data(kk,4)=sum(data(i1-n:i1-1,2))/n;
            f_data(kk,5)=sum(data(i1-n:i1-1,3))/n;
            f_data(kk,6)=sum(data(i1-n:i1-1,4))/n; 
           % f_data(kk,7)=sum(data(i1-n:i1-1,5))/n; 
            out_zero=0;
        end
    end
    
% find mean value every 1m interval    
n=0;
i=i1-1;
out_zero=1;
    while out_zero == 1 
        i=i+1;
        if data(i,1)>0.5+k-1 && data(i,1)<=0.5+k
          n=n+1;
              else if data(i,1)> 0.5+k 
                  k=k+1;
                  kk=kk+1;
                  f_data(kk,1)=time;
                  f_data(kk,2)=time;
                  f_data(kk,3)=-k+1;
                  f_data(kk,4)=sum(data(i-n:i-1,2))/n;
                  f_data(kk,5)=sum(data(i-n:i-1,3))/n;
                  f_data(kk,6)=sum(data(i-n:i-1,4))/n;
                 % f_data(kk,7)=sum(data(i-n:i-1,5))/n;
                  
                  n=0;
                  i=i-1;
              end
        end
        
        if i==(leng-1);
            
            k=k+1;
            kk=kk+1;
            
            f_data(kk,1)=time;
            f_data(kk,2)=time;
            f_data(kk,3)=-k+1;
            f_data(kk,4)=sum(data(i-n+1:i+1,2))/(n+1);
            f_data(kk,5)=sum(data(i-n+1:i+1,3))/(n+1);
            f_data(kk,6)=sum(data(i-n+1:i+1,4))/(n+1);
           % f_data(kk,7)=sum(data(i-n+1:i+1,5))/(n+1);
            out_zero=0;
            break
        end
    end
end

%calculating density-------------------------------------------------------

tem=f_data(:,4);sal=f_data(:,6);
f_data(:,7)=999.842594+6.793952.*10.^(-2).*tem-9.095290*10^(-3).*tem.^(2)...
+1.001685*10^(-4).*tem.^(3)-1.120083*10^(-6).*tem.^(4)+6.536332*10^(-9).*tem.^(5)...
+8.24493*10^(-1).*sal-4.0899*10^(-3).*tem.*sal+7.6438*10^(-5).*tem.^2.*sal...
-8.2467*10^(-7).*tem.^3.*sal+5.3875*10^(-9).*tem.^(4).*sal-5.72466*10^(-3).*sal.^(3/2)...
+1.0227*10^(-4).*tem.*sal.^(3/2)-1.6546*10^(-6).*tem.^(2).*sal.^(3/2)+4.8314*10^(-4).*sal.^2; 
%--------------------------------------------------------------------------
    
leng=length(f_data);
td=2;tem_dist=f_data(1,1:3);
for i=1:1:leng  
    if f_data(i,3) == 0
       tem_dist(td,1:3)=f_data(i,1:3);
       td=td+1;
    end
    tem_dist(td-1,3)=spheric_dist(tem_dist(td-2,1),tem_dist(td-1,1),tem_dist(td-2,2),tem_dist(td-1,2));
    
    f_data(i,8)=sum(tem_dist(:,3));
    f_data(i,8)=f_data(i,8)/1000;
end

[f, p]=uiputfile('*.dat','save line data');
sfiledir=[p,f];
fid=fopen(sfiledir,'w');
fprintf(fid,'%6d %6d %3d %8.4f %8.4f  %8.4f %8.4f %8.4f\r\n',f_data');
fclose(fid);

