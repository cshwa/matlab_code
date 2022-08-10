



close all
clear all
clc

%%%%%%%%%%%%%%%%%% load position OSTIA (1/20 degree) %%%%%%%%%%%%%%%%%%%%%%
% np=netcdf('E:\Downloads\ostia\2006\20060101-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA.nc');
% lona=np{'lon'}(:);
% lata=np{'lat'}(:);
% close(np)
% st_loa = find(lona <= 123.03 & lona >= 123);
% ed_loa = find(lona <= 133.03 & lona >= 133);
% st_laa = find(lata <=  30.03 & lata >=  30);
% ed_laa = find(lata <=  40.03 & lata >=  40);
% lona=lona(st_loa:ed_loa);
% lata=lata(st_laa:ed_laa);
% [lona lata] = meshgrid(lona,lata);

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

path_avhrr = 'f:\avhrr\'; %E:\avhrr\2010\AVHRR\
head_avhrr='avhrr-only-v2.';
% path_ostia = 'E:\Downloads\ostia\'; %E:\Downloads\ostia\2010\
% head_ostia='-UKMO-L4HRfnd-GLOB-v01-fv02-OSTIA';
foot='.nc';

st_year=1984;
ed_year=2013;

for year=st_year:ed_year
    ye=num2char(year,4);
    dd=yeardays(year);
    if dd==366
        days = [31 29 31 30 31 30 31 31 30 31 30 31];
    else
        days = [31 28 31 30 31 30 31 31 30 31 30 31];
    end
%   
%     month = 2:2:12;
%     m_len=length(month);

    for mon=2:2:12
        mo=num2char(mon,2);
        filek=(['.\re',ye,'\',ye,mo,'.dat']);
        fileout=(['.\re',ye,'\',ye,mo,'_comp.txt']);
        fid=fopen(fileout,'w+');
        kodc=load(filek);
        [kx ky] = size(kodc);
        stkk=kodc(:,1);
        yekk=kodc(:,2);
        monk=kodc(:,3);
        dayk=kodc(:,4);
        lonk=kodc(:,5);
        latk=kodc(:,6);
        temk=kodc(:,8);
%         for day=1:days(mon)
%             da=num2char(day,2);
            
%             % OSTIA file load
%             filea=([path_ostia,ye,'\',ye,mo,da,head_ostia,foot]);
%             nca=netcdf(filea);
% 
%             ssta=nca{'analysed_sst'}(:);
%             ssta=squeeze(ssta(st_laa:ed_laa,st_loa:ed_loa));
%             ssta(find(ssta==-32768))=NaN;
%             ssta=ssta*0.00999999977648258;
% %             maska=nca{'mask'}(:);
%             close(nca);
% %             maska=squeeze(maska(st_laa:ed_laa,st_loa:ed_loa));
% %             warning off
% %             maska=maska./maska;
% %             warning on
% %             ssta=ssta.*maska;
%             tema=griddata(lona,lata,ssta,lonk,latk);
%             tema_nan = find(isnan(tema));
%             tema(tema_nan) = -10000;
%             difa=temk-tema;
%             
            % avhrr file load
%             np=netcdf('E:\avhrr\2010\AVHRR\amsr-avhrr-v2.20100101.nc');
%             if mon == 12
%                 for jj = 1:2
%                     aa = year+jj-1;
%                     aa = year;
                    
                
            
                for kk = 1:kx
                    if yekk(kk) >= 2003 && yekk(kk)<= 2010
                        head_avhrr='amsr-avhrr-v2.';
                    else
                        head_avhrr='avhrr-only-v2.';
                    end
                    filer=([path_avhrr,num2char(yekk(kk),4),'\AVHRR\',head_avhrr,num2char(yekk(kk),4),num2char(monk(kk),2),num2char(dayk(kk),2),foot]);
                    ncr=netcdf(filer);
            
                    sstr=ncr{'sst'}(:);
                    close(ncr);
                    sstr=squeeze(sstr(st_lar:ed_lar,st_lor:ed_lor));
                    sstr(find(sstr==-999))=NaN;
                    sstr=sstr*0.00999999977648258;
                    temr=griddata(lonr,latr,sstr,lonk,latk);
                    temr_nan = find(isnan(temr));
                    temr(temr_nan) = -10000;
                    difr=temk-temr;
%                     if aa == yekk(kk)
%                         if day == dayk(kk)
                            if temr(kk) >= -100 && temr(kk) >= -100
                                disp([yekk(kk),monk(kk),dayk(kk)])
                                fprintf(fid,'%5d %5d %3d %3d %12.7f %12.7f %10.3f %10.3f %10.3f\n',stkk(kk),yekk(kk),monk(kk),dayk(kk),lonk(kk),latk(kk),temk(kk),temr(kk),difr(kk));
                            end
%                         end
%                     end
                end
%                 end
% %             else
%                 for kk = 1:kx
%                     aa=year;
%                     if aa >= 2003 && aa <= 2010
%                       head_avhrr='amsr-avhrr-v2.';
%                     else
%                       head_avhrr='avhrr-only-v2.';
%                     end
%                     filer=([path_avhrr,num2char(yekk(kk),4),'\AVHRR\',head_avhrr,num2char(aa,4),num2char(monk(kk),2),num2char(dayk(kk),2),foot]);
%                     ncr=netcdf(filer);
%             
%                      sstr=ncr{'sst'}(:);
%                      close(ncr);
%                     sstr=squeeze(sstr(st_lar:ed_lar,st_lor:ed_lor));
%                     sstr(find(sstr==-999))=NaN;
%                     sstr=sstr*0.00999999977648258;
%                     temr=griddata(lonr,latr,sstr,lonk,latk);
%                     temr_nan = find(isnan(temr));
%                     temr(temr_nan) = -10000;
%                     difr=temk-temr;
%                     if aa == yekk(kk)
% %                         if day == dayk(kk)
%                             if temr(kk) >= -100 && temr(kk) >= -100
%                                 disp([yekk(kk),monk(kk),dayk(kk)])
%                                 fprintf(fid,'%5d %3d %3d %12.7f %12.7f %10.3f %10.3f %10.3f\n',stkk(kk),monk(kk),dayk(kk),lonk(kk),latk(kk),temk(kk),temr(kk),difr(kk));
%                             end
% %                         end
%                     end
%                 end
%             end
            
        end
        fclose(fid);
end

% end
            
            
            
