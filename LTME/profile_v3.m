clear all;close all;clc
filepath = 'D:\장기생태\관측\20190823\CTD_Idronaut\alter\';
% cast = 2 %1=down , 2=up;

% ccf = [10 14 11 13 12 5 4 3 2 1 8 9 6 7]; %% snu point num.
% line =[7 14 8 12 9 3 2 4 13 10 5 11 1 6];  %RIST num
ccf = [10 14 11 13 12 5 4 3 2 1 8 9 6 7]; %% snu point num.
[out indx]= sort(ccf); % sort on snu point num.
line =[7 14 8 12 9 3 2 4 13 10 5 11 1 6];  %RIST num
line_sort = line(indx); % sort on Rist num
%%-------------------------------------------------------------------------
% % out = [1,2,3,4,5,6,7,8,9,10,11,12,13,14;]; %% snu point num.
% % line_sort == [10,13,4,2,3,1,6,5,11,7,8,9,12,14;]; %% RIST num
%%-------------------------------------------------------------------------
            
for i=1:14
    clearvars -except i out line_sort filepath  
    ki=line_sort(i);
    filename=[filepath,'CAST',num2str(i,'%04d'),'.txt'];
    fid=fopen(filename);
    if fid < 0
        continue
    else
    c_header=textscan(fid,'%s %s %s %s %s %s',1);
    c_data=textscan(fid,'%s %s %f %f %f %f');
    data=[c_data{3},c_data{4},c_data{6}];
    maxdepth=max(data(:,1));
    check_dep=find(maxdepth == data(:,1));
    data=data(1:check_dep,:);
    data=sortrows(data,1);
    leng=length(data);
    kk=0;
    k=0;
    n=0;
    out_zero=1;
    i1=0;
    while out_zero == 1
     i1=i1+1;
        if data(i1,1) <= 0.5
            n=n+1;
        else
            k=k+1;
            kk=kk+1;
            f_data(kk,1)=-floor(data(i1,1));
            f_data(kk,2)=sum(data(i1-n:i1-1,2))/n;
            f_data(kk,3)=sum(data(i1-n:i1-1,3))/n; 
            out_zero=0;
        end
    end
    
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
                  f_data(kk,1)=-k+1;
                  f_data(kk,2)=sum(data(i-n:i-1,2))/n;
                  f_data(kk,3)=sum(data(i-n:i-1,3))/n;
                  n=0;
                  i=i-1;
              end
        end
        
        if i==(leng-1);
             k=k+1;
             kk=kk+1;
            f_data(kk,1)=-k+1;
            f_data(kk,2)=sum(data(i-n+1:i+1,2))/(n+1);
            f_data(kk,3)=sum(data(i-n+1:i+1,3))/(n+1);
            out_zero=0;
            break
        end
    end
    end
    depth=f_data(:,1)*-1;
    temp=f_data(:,2);
    salt=f_data(:,3);
    dens=[];
%     for ii=1:1:length(depth)
%         dens(ii,1)=EOS80(temp(ii,1),salt(ii,1),depth(ii,1)*10);
%     end
%calculating density-------------------------------------------------------
% 
tem=temp;sal=salt;
dens=999.842594+6.793952.*10.^(-2).*tem-9.095290*10^(-3).*tem.^(2)...
+1.001685*10^(-4).*tem.^(3)-1.120083*10^(-6).*tem.^(4)+6.536332*10^(-9).*tem.^(5)...
+8.24493*10^(-1).*sal-4.0899*10^(-3).*tem.*sal+7.6438*10^(-5).*tem.^2.*sal...
-8.2467*10^(-7).*tem.^3.*sal+5.3875*10^(-9).*tem.^(4).*sal-5.72466*10^(-3).*sal.^(3/2)...
+1.0227*10^(-4).*tem.*sal.^(3/2)-1.6546*10^(-6).*tem.^(2).*sal.^(3/2)+4.8314*10^(-4).*sal.^2; 
%--------------------------------------------------------------------------
    depth=-depth;
    xlabels{1} = 'Temperature (^{˚}C)';
    xlabels{2} = 'Salinity(psu)';
    ylabels{1} = 'Depth(m)';
    ylabels{2} = 'Depth(m)';
    figure('position',[400 100 501 551],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.9]);
    set(gca,'Position',[0.10 0.06 0.75 0.85]);
    hold on;
    [ax,hl1,hl2]=plotxx(temp,depth,salt,depth,xlabels,ylabels);
    title_name=['Station-',num2str(ki,2)];
    title(title_name,'fontsize',17)
    set(hl1,'linewidth',1.5);set(hl2,'linewidth',1.5);
    set(get(ax(1),'ylabel'),'fontsize',16)
    set(get(ax(1),'xlabel'),'fontsize',16)
    set(get(ax(2),'ylabel'),'fontsize',16)
    set(get(ax(2),'xlabel'),'fontsize',16)
    set(ax(2),'fontsize',15,'xlim',[29 33.5],'ylim',[-30 0])
    set(ax(1),'fontsize',15,'xlim',[5 35],'ylim',[-30 0])
    out_name=['plot_v3_',num2str(ki,3)];
    saveas(gcf,out_name,'tif');
    
    figure('position',[400 100 501 551],'PaperUnits','inches','PaperPosition',[0 0 5.7 5.9]);
    set(gca,'Position',[0.15 0.10 0.70 0.80]);
    hold on
    title_name=['Station-',num2str(ki,2)];
    title(title_name,'fontsize',17)
    plot(dens-1000,depth,'g','linewidth',1.5);
    ylim([-30 0]);xlim([15 28]);
    xlabel('Density(σ)','fontsize',16)
    ylabel('Depth(m)','fontsize',16)
    set(gca,'box','on','fontsize',16)
   out_name=['plot_d_v3_',num2str(ki,3)];
    saveas(gcf,out_name,'tif');
    close all
    end



