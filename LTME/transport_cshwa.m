clear all; close all;
filename = 'J:\장기생태_2021\Dynamic\result\3regime\ncra_2reg_cir_2yr.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('GY_cshwa');
selectdepth=500;

s_max = 35.1000;
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
%% load picked points 
load picked_transport_cal_points.mat

t_point= 6;  %test1
% point_name={'r_mouth(low)', 'r_mouth(high)','l_mouth', 'myodo(high)', 'myodo(low)', 'left(ent)', 'main_channel', 'jinju_channel', 'small_posco_channel'};
point_name={'r_mouth(low)', 'r_mouth(high)','l_mouth', 'left(ent)', 'main_channel', 'jinju_channel'};
% kts = [p1_1 p1_2 ...  
%                p2_1 p2_2  ...
%                p3_1 p3_2  ...
%                p4_1 p4_2  ...
%                p5_1 p5_2  ...
%                p6_1 p6_2  ...
%                p7_1 p7_2  ...
%                p8_1 p8_2  ...
%                p9_1 p9_2];

kts = [p1_1 p1_2 ...  
               p2_1 p2_2  ...
               p3_1 p3_2  ...
               p6_1 p6_2  ...
               p7_1 p7_2  ...
               p8_1 p8_2];
    %% cshwa from editmask


depth = 100;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_1.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%r_mouth(low) r_mouth(high) l_mouth left(ent) main_channel jinju_channel\n');
    fprintf(fid,'%12.10f %12.10f %12.10f %12.10f %12.10f %12.10f %12.10f\n', trans');
%     fprintf(fid,'%%r_mouth(low) r_mouth(high) l_mouth myodo(high) myodo(low) left(ent) main_channel jinju_channel small_posco_channel\n');
%     fprintf(fid,'%8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f\n', trans');
    fclose(fid);
    
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_14.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];


% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_2.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_15.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];


% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_3.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_16.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];


% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_4.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_17.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];


% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];
% 

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_5.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_18.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};
kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];


% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_6.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_19.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];


% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_7.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_20.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];


% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_8.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_21.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];
% 

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_9.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_22.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];
% 

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_10.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_23.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];
% 

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_11.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all;
filename = 'monthly_spinup_24.nc';

clc; tic;warning off
foot = '.nc';
yyyy=2007;
year_start = yyyy;
year_end   = yyyy;
mon_start  = 1;
mon_end    =12;

g = grd('NWP');
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

t_point= 3;  %test1
point_name={'korea','tsugaru','soya', 'yellowsea'};

kts = [128.5000 34.8220 129.5000 33.1631 ...  
               140.4000  41.7249  140.3000  40.9738  ...
    141.9000  45.9794 141.8000  45.2092];

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];
% 

% kts = [128.883481062 35.2833167966 130.094555939 33.4263348468 ...  
%                140.1045  41.6281  140.2559  40.6754  ...
%     142.151473264  46.4790236999 142.043821201  45.1872096987];



depth = 1000;
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
    hz=g.z_w(2:N+1,:,:)-g.z_w(1:N,:,:); % z_w: [31x142x254]
    dx = 1./g.pm;
    dy = 1./g.pn;
    dx_v=0.5*(dx(1:M-1,:)+dx(2:M,:));
    dy_u=0.5*(dy(:,1:L-1)+dy(:,2:L));
    g_h_v=0.5*(g.h(1:M-1,:)+g.h(2:M,:));
    g_h_u=0.5*(g.h(:,1:L-1)+g.h(:,2:L));

    avgdxdy = mean([ mean( mean( dx ) ), mean( mean( dy ) ) ]);

    ynum = 0; mnum=0;
%     for yy = year_start : year_end ;
%         org_file=['..\',num2str(yy),'\monthly_mean\ocean_avg_mon_'];
%         org_file=['..\monthly\spinup3\ocean_avg_mon_'];
%         y_num = num2str(yy); 
            ynum=ynum+1;
%         for mm =  mon_start : mon_end
            mnum=mnum+1;
%             if mm<10
%                 mid=['0',num2str(mm)];
%             else
%                 mid=[num2str(mm)];
%             end

            nc=netcdf(filename,'read');
            zeta=nc{'zeta'}(:); % zeta(time, eta_rho, xi_rho)
            u=nc{'u'}(:);       % u(time, s_rho, eta_u, xi_u)
            v=nc{'v'}(:);       % v(time, s_rho, eta_u, xi_u)
            temp=nc{'temp'}(:); % temp(time, s_rho, eta_rho, xi_rho)
            %salt=nc{'salt'}(:); % salt(time, s_rho, eta_rho, xi_rho)
            close(nc)
          
            u(find(u<-1000))=0;
            v(find(v<-1000))=0;
            zeta(find(zeta<-1000))=0;
            temp(find(temp<-1000))=0;
            
            u = u.*mask3d_u;
            v = v.*mask3d_v;
            temp = temp.*mask3d_rho; % zero for land, ** very important **
%             salt = salt.*mask3d_rho; % zero for land, ** very important **
           
            %szero = ones(size(temp)).*s_max;
            %fresh = ( szero - salt ) ./ szero; % freshwater fraction
            %fresh = fresh.*mask3d_rho; % zero for land, ** very important **

            %   vertical coordinate changes in time
            %   because sea surface height changes in time.
            %   thickness of each layer changes propotional to total water thicknes.

            h_total = g.h + zeta;       %total water thickness
            for level=1:N               %thickness of each layer
                Hz(level,:,:)=squeeze(hz(level,:,:)).*(h_total./g.h);
            end

            % average Hz to  Arakawa-C u points

            Hz_u=0.5*(Hz(:,:,1:L-1)+Hz(:,:,2:L)); % each layer thickness
            z_u(1,:,:)=-g_h_u(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_u(k,:,:)=z_u(k-1,:,:)+Hz_u(k-1,:,:);
            end
            
            temp_u=0.5*(temp(:,:,1:L-1)+temp(:,:,2:L)); % each layer temp at u point
            %salt_u=0.5*(salt(:,:,1:L-1)+salt(:,:,2:L)); % each layer salt at u point
            %fresh_u=0.5*(fresh(:,:,1:L-1)+fresh(:,:,2:L)); % each layer freshwater at u point
            
                % average Hz to  Arakawa-C v points

            Hz_v=0.5*(Hz(:,1:M-1,:)+Hz(:,2:M,:)); % each layer thickness
            z_v(1,:,:)=-g_h_v(:,:);             % z @ bottom of each layer
            for k=2:+1:N
                z_v(k,:,:)=z_v(k-1,:,:)+Hz_v(k-1,:,:);
            end
                temp_v=0.5*(temp(:,1:M-1,:)+temp(:,2:M,:)); % each layer temp at u point
                %salt_v=0.5*(salt(:,1:M-1,:)+salt(:,2:M,:)); % each layer salt at u point
                %fresh_v=0.5*(fresh(:,1:M-1,:)+fresh(:,2:M,:)); % each layer freshwater at u point

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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=0;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
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
                            plot(  xpoint,  ypoint , 'ro')
                            dist2endpoint(icount) = m_lldist([xpoint xone],[ypoint yone])*1000;
                        else
                            tmp_dist=dist(icount);
                            vpoint=1;
                            icount=icount-1;
                            plot(  xpoint,  ypoint , 'gx')
                        end

                    end % if ( on_vpoint == 1 )

                    if( icount > 3 &&  dist2endpoint(icount) > dist2endpoint(icount-3) )
                        flag_approach = 0;
                    end
                end % while

            end % if (west2east_transect)

    total_temp=0;
    total_volume=0;
%     total_salt=0;
%     total_fresh=0;

            for index=1:icount

                vpoint=on_vpoint(index);
% % %                 calculate transport using xy_transport_function
                xy_transport_function

                if( west2east_transect == 0 && vpoint == 1 )
                    total_volume = total_volume - sum_segment;
%                     total_salt   = total_salt   - sum_segment_salt;
                    total_temp  = total_temp  - sum_segment_temp;
%                     total_fresh  = total_fresh  - sum_segment_fresh;
%                     freshtransect(index)=-sum_segment_fresh;
%                     voltransect(index)=-sum_segment;
                else
                    total_volume = total_volume + sum_segment;
%                     total_salt   = total_salt   + sum_segment_salt;
                    total_temp  = total_temp  + sum_segment_temp;
%                     total_fresh  = total_fresh  + sum_segment_fresh;
%                     freshtransect(index)=sum_segment_fresh;
%                     voltransect(index)=sum_segment;
                end


            end
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),' transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
%            disp(['YY/MM = ',num2str(yy),'/',num2char(mm,2),'heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])
           disp(['transport ',char(point_name(st)),': ',num2str(total_volume/1e+6)])
           disp(['heat transport ',char(point_name(st)),': ',num2str(total_temp*(4.1*10^6)/1e+15)])

           
           trans(mnum,st)    = total_volume/1e+6; 
            %             salt_tr(mnum,st)  = total_salt*1.025/1e+6;
                        temp_tr(mnum,st)  = total_temp*(4.1*10^6)/1e+15;
%         end   % for mm
%     end %% for yy
end  % for st=1:t_point

% salt = salt_tr ./ trans ;
%     trans(:,7)=trans(:,7)*(-1);
%     outfile = ['roms_transport_monthly.txt'];
    outfile = ['roms_transport_monthly_csh_12.txt'];
    fid = fopen(outfile,'w+');
    fprintf(fid,'%%korea tsugaru soya\n');
    fprintf(fid,'%8.3f %8.3f %8.3f\n', trans');
    fclose(fid);
clear all; close all; clc;
