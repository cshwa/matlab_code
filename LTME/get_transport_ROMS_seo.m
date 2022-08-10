% function [total_volume,total_salt,total_fresh] = xy_transport_onefile(filename,g,selectdepth);
% ====================================================================
% [total_volume,total_salt,total_fresh]= xy_transport_onefile(filename,g,selectdepth);
% calculate transports of  volume(in Sv, 10^6 m^3/sec)
%                          heat  (in PW, 10^15 W)
%                          salt  (in 10^9 kg)
%
% across a line segment (slice along a constant I or J)
%
% 'filename'       = history or average file name
% 'grid_nick_name' = grid information such as 'eas' 'hudson' 'latte'
% 'selectdepth' = surface to which depth (m) such as 100, 500, 50000
%
% keep in mind that
% vertical coordinate changes in time = h + zeta(t) in ROMS
%
% Updated 14-Apr-2018 by Yong-Yub Kim.
% ====================================================================

clc;close all;clear all;
warning off;

linux=1; windows=0;

if (windows ==1)
    % % for windows
    dropboxpath='C:\Users\KYY\Dropbox'; %% SNU_desktop
    addpath(genpath([dropboxpath '\source\matlab\Common\m_map']));
    addpath(genpath([dropboxpath '\source\matlab\Common\Figure']));
    addpath(genpath([dropboxpath '\source\matlab\Common\netcdf_old']));
    addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Grid_kyy']));
%     addpath(genpath([dropboxpath '\source\matlab\Model\ROMS\Analysis\Figure\seo_nwp_1_10\run']));
elseif (linux==1)
    % % for linux
%     dropboxpath='/home01/kimyy/Dropbox'; %% ROMS server
    dropboxpath='/home/kimyy/Dropbox'; %% DAMO server
    addpath(genpath([dropboxpath '/source/matlab/Common/m_map']));
    addpath(genpath([dropboxpath '/source/matlab/Common/Figure']));
    addpath(genpath([dropboxpath '/source/matlab/Common/netcdf_old']));
    addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Grid_kyy']));
%     addpath(genpath([dropboxpath '/source/matlab/Model/ROMS/Analysis/Figure/seo_nwp_1_10/run']));
end



% % for dgist(cannot use dropbox code in the dgist)
% % addpath(genpath(['/scratch/snu01/kimyy/roms_nwp/nwp_1_20/output/',testname,'/TRANSP/m_map']))
% % addpath(genpath(['/scratch/snu01/kimyy/roms_nwp/nwp_1_20/output/',testname,'/TRANSP/netcdf_old']))


yearstart = 1980;
yearend = 2009;
testname = 'avg_ens_10km_mean_monthly';
foot = '.nc';
inputdir = '/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_monthly/'; %% for DAMO, nwp_seo_1/10
outputdir = '/data2/kimyy/Reanalysis/nwp_1_10_seo/ens_mean_monthly/TRANSP/';  %% for DAMO, nwp_seo_1/10
outfile = strcat(outputdir,'roms_tr_nwp_seo_monthly_',num2str(yearstart),'_',num2str(yearend),'_grid_10km_new.txt');

fid = fopen(outfile,'w+');
fprintf(fid,'%%korea tsugaru soya tartar\n');
fclose(fid);

for year = yearstart:yearend
    for month=1:12
%     prefix = '/scratch/snu01/kimyy/roms_nwp/nwp_1_20/output/';
    %prefix = '/scratch/snu02/roms_nwp/output/smooth13_vtvs/yearly_'
    tempyear = num2str(year,'%04i');
    tempmonth = num2str(month,'%02i');
    %numb ='08';
    filename = strcat(inputdir,testname,'_',tempyear,'_',tempmonth,foot)

%     yyyy=2001;
%     year_start = yyyy;
%     year_end   = yyyy;
    % mon_start  = 1;
    % mon_end    =12;
%     if (exist('g' , 'var') ~= 1)
        g = grd('seo_NWP_reanalysis');
%     end
    selectdepth=5000;

    s_max = 34.573;
    % *************************************************************
    %
    %   END OF USER DEFINED VARIABLES
    %
    % *************************************************************

    % size of grids
    [r,c] = size ( g.lon_rho );
    mask3d_rho=repmat(g.mask_rho,[1 1 g.N]);
    mask3d_rho=permute(mask3d_rho,[3 1 2]);

    [r1,c1] = size ( g.lon_u );
    mask3d_u=repmat(g.mask_u,[1 1 g.N]);
    mask3d_u=permute(mask3d_u,[3 1 2]);

    [r2,c2] = size ( g.lon_v );
    mask3d_v=repmat(g.mask_v,[1 1 g.N]);
    mask3d_v=permute(mask3d_v,[3 1 2]);

    % transport from surface to which depth (m)

    if ( selectdepth > 0 )
        selectdepth = selectdepth*-1;
    end

    t_point= 4;  %test1
    point_name={'korea','tsugaru','soya', 'tartar'};

    kts = [128.0 35.2 134.0 35.2 ...  
           140.1  41.6  140.3  40.5  ...
        142.2  46.4 142.0  45.0 ...
        138.5  47.5 142.2  47.5
        ];

    depth = 5000;
    sn=0;

    for st=1:t_point

        disp([char(point_name(st)),' Depth to ',num2str(depth)])
        endpt_lon(1) = kts(sn+1);
        endpt_lat(1) = kts(sn+2);
        endpt_lon(2) = kts(sn+3);
        endpt_lat(2) = kts(sn+4);
        sn=sn+4;

        for cpts=1:2 % corner points
    % %         find points that are close with stations 
            dist = sqrt(  ( g.lon_rho - endpt_lon(cpts) ).*( g.lon_rho - endpt_lon(cpts) ) + ...
                ( g.lat_rho - endpt_lat(cpts) ).*( g.lat_rho - endpt_lat(cpts) ) );
            ind=find( min( dist(:) ) == dist );
            % closest points row and column indice
            row_index = mod ( ind - 1, r ) + 1;
            col_index = floor( (ind - 1) / r ) + 1;
            corner_endpt_col(cpts)=col_index(1);
            corner_endpt_row(cpts)=row_index(1);
        end


        % my xy_transport_onefile works only if corner_endpt_row(2) < corner_endpt_row(1).
        if( corner_endpt_row(2) > corner_endpt_row(1)  )
            tmp_col=corner_endpt_col(2);
            tmp_row=corner_endpt_row(2);
            corner_endpt_col(2)=corner_endpt_col(1);
            corner_endpt_row(2)=corner_endpt_row(1);
            corner_endpt_col(1)=tmp_col;
            corner_endpt_row(1)=tmp_row;
            beep
            disp(' === switching two end points === ')
        end

        % longitude and latitude coordinate.

        for i=1:length(endpt_lat)
            xx(i)=g.lon_rho(corner_endpt_row(i),corner_endpt_col(i));
            yy(i)=g.lat_rho(corner_endpt_row(i),corner_endpt_col(i));
        end
        distance_r = m_lldist ( xx, yy )*1000;

        %  transect information

        if( corner_endpt_col(2) >= corner_endpt_col(1) )
            delj=1;
            west2east_transect=1; % previously zonaltransect
        else
            delj=-1;
            west2east_transect=0; % previously meridionaltransect
        end

        if( corner_endpt_row(2) > corner_endpt_row(1) )
            deli=1;
        else
            deli=-1;
        end

        xzero=g.lon_rho( corner_endpt_row(1), corner_endpt_col(1) );
        yzero=g.lat_rho( corner_endpt_row(1), corner_endpt_col(1) );
        xone=g.lon_rho( corner_endpt_row(2), corner_endpt_col(2) );
        yone=g.lat_rho( corner_endpt_row(2), corner_endpt_col(2) );
        slope=( yone-yzero) / (xone - xzero);
        % A x + B y + C = 0;
        A=slope;
        B=-1;
        C=-slope*xzero+yzero;
        D=sqrt( A*A + B*B );
        % distance = abs( A x + B y + C ) / D

        %   grid information

        %N  is the number of vertical levels
        %hz is thickness  of each level
        N = g.N;
        [M L]=size(g.h);
        if (g.Vtransform==1)
            hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
        elseif(g.Vtransform==2)
% %         hz won't be used in Vtransform =2
    %        for i=1:N
    %            hz(i,:,:)=squeeze(g.z_w(i+1,:,:))./g.h - squeeze(g.z_w(i,:,:))./g.h;  %% Y.Y.KIM            
    %        end
        end
        dx = 1./g.pm;
        dy = 1./g.pn;
        dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
        dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
        g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
        g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

        avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

        ynum = 0; mnum=0;
                ynum=ynum+1;
                mnum=mnum+1;

                nc=netcdf(filename,'read');
                zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
                u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
                v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
                temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
                close(nc)
% 
                u(find(u<-1000))=0;
                v(find(v<-1000))=0;
                zeta(find(zeta<-1000))=0;
                temp(find(temp<-1000))=0;

                u = u.*mask3d_u;
                v = v.*mask3d_v;
                temp = temp.*mask3d_rho; % zero for land, ** very important **

                %   vertical coordinate changes in time
                %   because sea surface height changes in time.
                %   thickness of each layer changes propotional to total water thicknes.

                if (g.Vtransform==1)  %%Y.Y.KIM.
                  h_total = g.h + zeta;       %total water thickness
                elseif(g.Vtransform==2)
                    for i=1:N+1
                        S_w_nossh(i,:,:)=squeeze(g.z_w(i,:,:))./g.h;
                        g.z_w(i,:,:)=zeta+(zeta+g.h).* squeeze(S_w_nossh(i,:,:));
                    end
                end
                for level=1:N               %thickness of each layer
                    if (g.Vtransform==1)            
                     %% hz -> layer thickness when zeta = 0 (grid file)
                     %% Hz -> zeta in the result file is reflected
                        Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
                    elseif(g.Vtransform==2)
                        Hz(level,:,:)=g.z_w(level+1,:,:)-g.z_w(level,:,:); %%%Y.Y.KIM.
                    end  
                end

                % average Hz to  Arakawa-C u points

                Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
                z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
                for k=2:+1:N
                    z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
                end

                temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point

                    % average Hz to  Arakawa-C v points

                Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
                z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
                for k=2:+1:N
                    z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
                end
                    temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point

                %   ====================================================================================
                %   find path from corner_endpt(1) to corner_endpt(2)

                icount=1;
                col_index(icount)=corner_endpt_col(1);
                row_index(icount)=corner_endpt_row(1);
                on_vpoint(icount)=0;
                vpoint=0;
                xpoint=g.lon_u( row_index(icount), col_index(icount) );
                ypoint=g.lat_u( row_index(icount), col_index(icount) );

                signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;
                tmp_dist=dist(icount);
                dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                flag_approach = 1;

                if (west2east_transect)

                    while ( dist2endpoint(icount) > avgdxdy  &&  flag_approach )

                        icount=icount+1;

                        if ( vpoint == 1 )

                            col_index(icount)=col_index(icount-1)+delj;
                            if ( on_vpoint(icount-1) == 1)
                                row_index(icount)=row_index(icount-1);
                            else
                                row_index(icount)=row_index(icount-1)+deli;
                            end
                            xpoint=g.lon_v( row_index(icount), col_index(icount) );
                            ypoint=g.lat_v( row_index(icount), col_index(icount) );
                            signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                            dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                            if ( signline(icount)*signline(icount-1) < 0  ...
                                    ||   dist(icount) <= dist(icount-1)       ...
                                    ||   dist(icount) <= tmp_dist                  )
                                tmp_dist=0;
                                on_vpoint(icount)=1;
    %                             plot(  xpoint,  ypoint , 'ro')
                                dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                            else
                                tmp_dist=dist(icount);
                                vpoint=0;
                                icount=icount-1;
    %                             plot(  xpoint,  ypoint , 'gx')
                            end

                        else % on upoint

                            col_index(icount)=col_index(icount-1);
                            if ( on_vpoint(icount-1) == 0)
                                row_index(icount)=row_index(icount-1)+deli;
                            else
                                row_index(icount)=row_index(icount-1);
                            end
                            xpoint=g.lon_u( row_index(icount), col_index(icount) );
                            ypoint=g.lat_u( row_index(icount), col_index(icount) );
                            signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                            dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                            if (      signline(icount)*signline(icount-1) < 0 ...
                                    ||   dist(icount) <= dist(icount-1)   ...
                                    ||   dist(icount) <= tmp_dist                       )
                                tmp_dist=0;
                                on_vpoint(icount)=0;
    %                             plot(  xpoint,  ypoint , 'ro')
                                dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                            else
                                tmp_dist=dist(icount);
                                vpoint=1;
                                icount=icount-1;
    %                             plot(  xpoint,  ypoint , 'gx')
                            end

                        end % if ( on_vpoint == 1 )

                        if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                            flag_approach = 0;
                        end

                    end % while

                else % if (west2east_transect)

                    while ( dist2endpoint(icount) > avgdxdy  &&  flag_approach )

                        icount=icount+1;

                        if ( vpoint == 1 )

                            if ( on_vpoint(icount-1) == 1)
                                col_index(icount)=col_index(icount-1)+delj;
                                row_index(icount)=row_index(icount-1);
                            else
                                col_index(icount)=col_index(icount-1);
                                row_index(icount)=row_index(icount-1)+deli;
                            end

                            xpoint=g.lon_v( row_index(icount), col_index(icount) );
                            ypoint=g.lat_v( row_index(icount), col_index(icount) );
                            signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                            dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                            if ( signline(icount)*signline(icount-1) < 0  ...
                                    ||   dist(icount) <= dist(icount-1)       ...
                                    ||   dist(icount) <= tmp_dist                  )
                                tmp_dist=0;
                                on_vpoint(icount)=1;
    %                             plot(  xpoint,  ypoint , 'ro')
                                dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                            else
                                tmp_dist=dist(icount);
                                vpoint=0;
                                icount=icount-1;
    %                             plot(  xpoint,  ypoint , 'gx')
                            end

                        else % on upoint

                            if ( on_vpoint(icount-1) == 0)
                                col_index(icount)=col_index(icount-1);
                                row_index(icount)=row_index(icount-1)+deli;
                            else
                                col_index(icount)=col_index(icount-1)+delj;
                                row_index(icount)=row_index(icount-1);
                            end


                            xpoint=g.lon_u( row_index(icount), col_index(icount) );
                            ypoint=g.lat_u( row_index(icount), col_index(icount) );
                            signline(icount)=(ypoint - yzero) - slope*(xpoint -xzero);
                            dist(icount)= abs( A*xpoint + B*ypoint + C ) / D;

                            if (      signline(icount)*signline(icount-1) < 0 ...
                                    ||   dist(icount) <= dist(icount-1)   ...
                                    ||   dist(icount) <= tmp_dist                       )
                                tmp_dist=0;
                                on_vpoint(icount)=0;
    %                             plot(  xpoint,  ypoint , 'ro')
                                dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                            else
                                tmp_dist=dist(icount);
                                vpoint=1;
                                icount=icount-1;
    %                             plot(  xpoint,  ypoint , 'gx')
                            end

                        end % if ( on_vpoint == 1 )

                        if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                            flag_approach = 0;
                        end
                    end % while

                end % if (west2east_transect)

        total_temp=0;
        total_volume=0;

                for index=1:icount

                    vpoint=on_vpoint(index);
    % % %                 calculate transport using xy_transport_function
                    xy_transport_function

                    if( west2east_transect == 0 && vpoint == 1 )
                        total_volume = total_volume - sum_segment;
                        total_temp  = total_temp  - sum_segment_temp;
                    else
                        total_volume = total_volume + sum_segment;
                        total_temp  = total_temp  + sum_segment_temp;
                    end


                end
               disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
               disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])


               trans(mnum,st)    = total_volume/1e+6; 
               temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
    end  
    
    fid = fopen(outfile,'a+');
    fprintf(fid,'%8.3f %8.3f %8.3f %8.3f \n', trans');
    fclose(fid);
    end
end
