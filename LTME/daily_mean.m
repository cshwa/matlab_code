



close all
clear all
clc

%%%%%%%%%%%%%%%%%% load position avhrr (1/4 degree) %%%%%%%%%%%%%%%%%%%%%%%%
np=netcdf('f:\avhrr\2006\AVHRR\amsr-avhrr-v2.20060101.nc');
lonr=np{'lon'}(:);
latr=np{'lat'}(:);
close(np)
st_lor = find(lonr <= 123.2 & lonr >= 123);
ed_lor = find(lonr <= 133.2 & lonr >= 133);
st_lar = find(latr <=  30.2 & latr >=  30);
ed_lar = find(latr <=  40.2 & latr >=  40);
lonr=lonr(st_lor:ed_lor);
latr=latr(st_lar:ed_lar);
[lonr latr] = meshgrid(lonr,latr);
[y, x] = size(lonr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = load('position_all2.dat');
lonk = data(:,2);
latk = data(:,3);
len_k = length(lonk);

path_avhrr = 'f:\avhrr\';
head_avhrr='avhrr-only-v2.';
foot='.nc';

st_year=1984;
ed_year=2013;

t_sst = zeros(365,len_k);
nn = zeros(365,1);
for year=st_year:ed_year
    ye=num2char(year,4);
    infile_name  = ([path_avhrr,ye,'_re.dat']);
    infile=textread(infile_name,'%25c');
    [was mm2 dd2]=textread(infile_name,'%18s%2d%2d');
    clear was;
    aa = size(mm2)/2;
    for i = 1:aa
        mm(i) = mm2(i*2-1);
        dd(i) = dd2(i*2-1);
    end
    clear mm2; clear dd2;
    
    days =  [31 28 31  30  31  30  31  31  30  31  30  31];
%     days =  [31 59 90 120 151 181 212 243 273 304 334 365];
    
    
    for i = 1:length(infile)
        infile2 = ([path_avhrr,num2char(year,4),'\AVHRR\',squeeze(infile(i,:))]); %'f:\avhrr\2006\AVHRR\
        ncr = netcdf(infile2);
        sstr=ncr{'sst'}(:);
        close(ncr);
        sstr=squeeze(sstr(st_lar:ed_lar,st_lor:ed_lor));
        sstr(find(sstr==-999))=NaN;
        sstr=sstr*0.00999999977648258;
        n = 0;
        for j = 1:12
            for k = 1:days(j)
                n = n+1;
                if j==mm(i)
                    if k==dd(i)
                        jd(i) = n;
                        nn(jd(i)) = nn(jd(i)) + 1;
                    end
                end
            end
        end
        sst = griddata(lonr,latr,sstr,lonk,latk);
        t_sst(jd(i),:) = t_sst(jd(i),:) + sst'; 
    end
    disp([num2char(year,4)])
end

fid = fopen('dailiy_mean.dat','w+');
for i = 1:365
    f_sst= nanmean(t_sst(i,:)/nn(i));
    qq = nn(i);
    fprintf(fid,'%10.3f%5d\n',f_sst, qq);
end
fclose(fid);
            
            
            
