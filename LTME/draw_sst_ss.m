% domaxis = [124 132 32 36];
domaxis=[126.125 129.125 33.125 35.125]; 
s_yr = 1982;
e_yr = 2016;
for yr = 2013
filepath1 = ['E:\11.사업\장기생태_3단계\1차년도\Data\위성자료\daily\' num2str(yr) '\'];%daily folder
filepath2 ='E:\11.사업\장기생태_3단계\1차년도\Data\위성자료\monthly\';   %monthly folder
name1 = ['avhrr-only-v2.' num2str(yr)];
name11 = ['avhrr_monthly' num2str(yr)];

    for i = 8
        name2=[num2char(i,2)];
%         if i < 10 
%             name2 = strcat('0', num2str(i));
%         else
%             name2 = num2str(i);
%         end
        name11 = ['avhrr_monthly' num2str(yr) '_' name2 '.nc'];
         %define name1(month name)

        file1 = strcat(filepath2, name11);
        nc = netcdf(file1);
        row_lat = nc{'lat'}(:); %unit: degrees_E  
        row_lon = nc{'long'}(:);%unit: degrees_N 
        row_temp= nc{'temp'}(:);
        clear nc;    

        figure('position',[100 100 500 400]);
        hold on  
        pcolor(row_lon,row_lat,row_temp);
    %     shading flat;
        shading interp
        title([num2str(yr),'. ', num2str(i)],'fontsize',14);
        caxis([10 30]);
        c=colorbar;
        set(gca, 'xlim',[domaxis(1) domaxis(2)],'ylim',[domaxis(3) domaxis(4)]);
    %     colorscale([0 64],[30 35],2,'vert','position',[0.92 0.15 0.03 0.7]);
        set(c,'YaxisLocation','right');
        ylabel(c,'SST (Deg.)') 

        %--- save figure --------------------
        outpath = 'figure\';
        out_vari = name11;    
%         out=[outpath,num2str(name2),'_',num2str(yr)];
        out=[outpath,'avhrr_monthly' num2str(yr) '_' name2];
        set(gcf,'renderer','painter');
        set(gcf, 'PaperUnits', 'inches');
        x_width=5*1.5;
        y_width=4*1.3;
        set(gcf, 'PaperPosition', [0 0 x_width y_width]); %
        saveas(gcf,out,'tif');
        %------------------------------------    
        clf;close all
    end        
end


