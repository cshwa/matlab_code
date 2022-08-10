
clear all; close all;

%==========================================================================
% gname= 'roms_grid_4degree_6.nc';
% gd       = netcdf(gname);
% lon_rho  = gd{'lon_rho'}(:);
% lat_rho  = gd{'lat_rho'}(:);
% mask_rho = gd{'mask_rho'}(:);
% close(gd)
% rmask = mask_rho./mask_rho;
% warning off;
%==========================================================================
% bln=load('coast_file.dat');
% korea_lat=bln(:,1);
% korea_lon=bln(:,2);
% % lim=[115 160 15 53]; % eas4
% lim=[117.5 154.5 18.4 48.5]; % eas10
%==========================================================================

in_path='..\yearly\'; 
mid_name = 'avhrr_'; foot ='.nc';
out_path = '..\SST_mean\vari\';
mkdir(out_path);

s_year=1982; e_year=2009;

merge_s=[]; merge2=[]; imerge=[];
ii=0;
for yid=s_year:e_year
     ii=ii+1;
 m_num = [num2str(yid)];
    filename = [in_path,mid_name,m_num,foot];
    nc=netcdf(filename);
    temp = nc{'temp'}(:);
    close(nc)

    tem_sur=temp.*rmask;
   
    bi=162*194;

    b_sur=reshape(tem_sur,bi,1);    
    merge_s=[merge_s b_sur];        
    end
%--------------------------------------------------------variation
res_s=[];  
for fi=1:bi
    x1=1:ii;  xxx=x1';
    p1_s=polyfit(x1,merge_s(fi,:),1);
    res_s=[res_s ; p1_s]; 
end

p2_s=res_s(:,1) ;
final_p_s=reshape(p2_s,162,194);

%%%%%%%%%%%%%%%%%%%%%%%%%%% surface plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
temp_lim = [-0.1 0.1];  level_a  = [-0.1:0.04:0.1];
color_c  = '-w';  out_type = 'tif';


pcolorjw(lon_rho,lat_rho,final_p_s);
hold on
[C,h]=contour(lon_rho,lat_rho,final_p_s,level_a,color_c,'linewidth',1.5);

clabel(C,h,'FontSize',10,'Color','k','labelspacing',100000,'Rotation',0,'fontweight','bold');
caxis(temp_lim );
geoshow(korea_lon, korea_lat,'DisplayType','polygon','facecolor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
geoshow(korea_lon, korea_lat,'DisplayType','polygon','facecolor',[0.7 0.7 0.7],'EdgeColor',[0.7 0.7 0.7]);
set(gca,'fontsize',18,'fontweight','bold','fontname','times new roman'); % X,Y축의 폰트
set(gca,'linewidth',2,'tickdir','in')
xlabel(['Longitude (^o E)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
ylabel(['Latitude(^o N)'],'fontsize',18,'fontweight','bold','fontname','times new roman');
title(['Satellite SST Trend '],'fontsize',18,'fontweight','bold');   %월바꾸기
        load('colorData.mat');
        set(gcf,'Colormap',Mtemperaturemap)
        axis(lim)
        hold off
%         colorbar
                colorscale([1 64],[-0.1, 0.1],0.04,'vert','position',[0.92 0.25 0.020 0.6]);
         set(gca,'YaxisLocation','right','fontsize',10,'fontname','times new roman');
out=[out_path,'sat_sst_trend'];
saveas(gcf,out,out_type)
close all

