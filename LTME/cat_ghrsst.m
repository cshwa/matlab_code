function C=cat_ghrsst(R)
% return structure C [Time,URLSITE,URLPATH,URLFILE]

    if isfield(R,'DataSetBranch') & strcmp(R.DataSetBranch,'GHRSST_SEVIRI_SST')
% catalog of files for GHRSST_SEVIRI_SST

% files have a regular name pattern with varying date and hour of observation
% example
%http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR/2008/092/
%'20090410-SEVIRI_SST-EUR-L2P-sst3mlml_20090410_0400-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_0400-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_0700-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1000-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1300-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1600-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_1900-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080402_2200-v01.nc.bz2';
%'20080402-SEVIRI_SST-EUR-L2P-sst3mlml_20080403_0100-v01.nc.bz2';
%
% permanent archive is at NODC, the latest few months at JPL

%J2=floor(now);J1=J2-30;
J1=floor(datenum(2005,01,31,13,00,00)); % beginning of dataset
J2=floor(now);                      
hs=[4 7 10 13 16 19 22 25]; % hours of several frames per day
Nt=(J2-J1+1)*length(hs); % number of frames = number of days * frames per day

% 1 frame per file
URLSITE=cell(Nt,1);
URLPATH=cell(Nt,1);
URLFILE=cell(Nt,1);
%URLCTIME=cell(Nt,1);
%TIME=cell(Nt,1);
Time=zeros(Nt,1);
%t0=datenum('1981-01-01 00:00:00','yyyy-mm-dd HH:MM:SS');

URLSITE1='http://data.nodc.noaa.gov/opendap/ghrsst/L2P/SEVIRI_SST/EUR'; % permanent archive
JURLSITE2=J2-45;
URLSITE2='http://dods.jpl.nasa.gov/opendap/sea_surface_temperature/GHRSST/data/L2P/SEVIRI_SST/EUR'; % most recent data (last few months)

kt=0;
        for j=J1:J2
        DATE1=datestr(j,'yyyymmdd');
        YYYY=datestr(j,'yyyy');
        yearday=(j-datenum(YYYY,'yyyy')+1);
        DDD=sprintf('%03d',yearday);
            for h=hs
            kt=kt+1;
            t=j+h/24;
Time(kt)=t;  
%TIME{kt}=datestr(t,'yyyy-mm-dd HH:MM:SS');
                if j<JURLSITE2
URLSITE{kt}=URLSITE1;
                else
URLSITE{kt}=URLSITE2; 
                end
URLPATH{kt}=['/' YYYY '/' DDD '/'];
            DATE2=datestr(t,'yyyymmdd_HHMM');        
URLFILE{kt}=[DATE1 '-SEVIRI_SST-EUR-L2P-sst3mlml_' DATE2 '-v01.nc.bz2'];
            end
        end
C.Time=Time;C.URLSITE=URLSITE;C.URLPATH=URLPATH;C.URLFILE=URLFILE;
    end %if R.DataSetBranch

