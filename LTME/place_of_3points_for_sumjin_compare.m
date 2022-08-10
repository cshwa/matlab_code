close all; clear all; clc;
cd D:\장기생태\Dynamic\06_river_physics\환경과학원
agyang=load('agyang_regression_yeosu_extract_climate_2021.mat'); 
hadong=load('hadong_regression_yeosu_extract_climate_2021.mat'); 
jinwal=load('jinwal_regression_yeosu_extract_climate_2021_high_airT.mat'); 
cd D:\장기생태\Dynamic\06_river\환경과학원
songjung=load('songjung_regression_yeosu_extract_climate_2021.mat');


% make 1980~present
clearvars ref_*
k=0
for i = 1989:2020
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:size(eom_d,1)
    l=0
        ref_yy(i,:)=[num2str(i+1988)];
    for n = 1:12
        m = m+1;
        ref_yymm{m}=[num2str(i+1988) '-' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        ref_yymmdd{k}=[num2str(i+1988) '-' num2str(n,'%02d') '-'  num2str(j,'%02d')];
        ref_mmdd{k}=[num2str(n,'%02d') '-'  num2str(j,'%02d')];
    end
    end
end

k=0; m=0;
for i = 1:size(eom_d,1)
    l=0
%         ref_yy(i,:)=[num2str(i+1988)];
    for n = 1:1
        m = m+1;
        ref_yy01{m}=[num2str(i+1988) '-' num2str(n,'%02d') '-'  num2str(1,'%02d')];
    end
end

% find matched date from yy mm dd form.
j=0
for i = 1:length(agyang.r_date_txt)
   if  sum(strcmp(agyang.r_date_txt(i,:), ref_yymmdd)) ~= 0
       j=j+1;
       indx_date_agyang(j) = find([strcmp(agyang.r_date_txt(i,:), ref_yymmdd)] == 1);    
   else
       disp(i)
   end
end

j=0
for i = 1:length(hadong.r_date_txt)
   if  sum(strcmp(hadong.r_date_txt(i,:), ref_yymmdd)) ~= 0
       j=j+1;
       indx_date_hadong(j) = find([strcmp(hadong.r_date_txt(i,:), ref_yymmdd)] == 1);    
   else
       disp(i)
   end
end

j=0
for i = 1:length(songjung.r_date_txt)
   if  sum(strcmp(songjung.r_date_txt(i,:), ref_yymmdd)) ~= 0
       j=j+1;
       indx_date_songjung(j) = find([strcmp(songjung.r_date_txt(i,:), ref_yymmdd)] == 1);    
   else
       disp(i)
   end
end

j=0
for i = 1:length(jinwal.r_date_txt)
   if  sum(strcmp(jinwal.r_date_txt(i,:), ref_yymmdd)) ~= 0
       j=j+1;
       indx_date_jinwal(j) = find([strcmp(jinwal.r_date_txt(i,:), ref_yymmdd)] == 1);    
   else
       disp(i)
   end
end

j=0
for i = 1:length(ref_yy01)
   if  sum(strcmp(ref_yy01{i}, ref_yymmdd)) ~= 0
       j=j+1;
       indx_date_x(j) = find([strcmp(ref_yy01{i}, ref_yymmdd)] == 1);    
   else
       disp(i)
   end
end


figure; hold on;
plot(indx_date_songjung,songjung.r_temp,'b','linew',2); 
plot(indx_date_hadong,hadong.r_temp,'g','linew',2); 
plot(indx_date_agyang,agyang.r_temp,'r','linew',2); 
plot(indx_date_jinwal,jinwal.r_temp,'k','linew',2); 
% plot(indx_date_agyang,agyang.r_temp,'ro'); 
% plot(indx_date_hadong,hadong.r_temp,'go'); 
% plot(indx_date_songjung,songjung.r_temp,'bo'); 
grid on;
xlim([indx_date_agyang(1) inf]);
xticks(indx_date_x);
xticklabels(1989:2020);
% legend('구례','하동','악양');
legend('구례','하동','악양','진월');
set(gca,'fontsize',16); ylabel('temperature (^oC)');
xlim([1 inf]);
xtickangle(45)




