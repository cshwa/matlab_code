%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clearvars;

var_name = { 'DO',   'Chla',   'NO3',   'NH4'};
dim_name = { 'mg/L', 'mg/m^3', 'mg/L', 'mg/L' };
var_col  = {   9,      16,       20,     21 };
flag_ud  = true;
st_name  = '구례_fix';
fig_name = '섬진(송정)-'
%%%%
% flag_ud  = true;
% st_name  = '사천천';
% fig_name = '사천천-'
%%%%
% flag_ud  = false;
% st_name  = '남강댐1';
% fig_name = '남강댐1-'
%%%%
% flag_ud  = true;
% st_name  = '하동';
% fig_name = '하동-'

x_tick = [];
j=0;
for i = 1989:1:2020
    j=j+1;
    x_tick = [x_tick, datenum( i,01,01 )] ;
     temp = datestr(datenum( i,01,01 ));
     x_tick_c{j,1} = temp(4:end);
end

cd E:\장기생태\Dynamic\06_river\환경과학원
[raw txt]=xlsread(['수질측정지점_',st_name,'.xls'],'검색결과','');

if( flag_ud ), r_txt_ud = flipud(txt);
else           r_txt_ud = txt;  end

%-- time 처리부분
if( flag_ud ), r_date_txt=r_txt_ud(1:end-1,5);
else           r_date_txt=r_txt_ud(2:end,5); end
time = zeros( length(r_date_txt), 1 ); 

for i = 1:length(r_date_txt)
    time(i) = datenum(r_date_txt{i},'yyyy.mm.dd');
    time_c{i} = datenum(r_date_txt{i},'yyyy.mm.dd');
end
%%
k=0
for i = 1989:2019
            k=k+1;
    for j = 1:12
        eom_d(k,j) = eomday(i,j); % 1980 is leap-yr
    end
end

k=0; m=0;
for i = 1:31
    l=0
        ref_yy(i,:)=[num2str(i+1988)];
    for n = 1:12
        m = m+1;
        ref_yymm(m,:)=[num2str(i+1988) '.' num2str(n,'%02d')];
    for j = 1:eom_d(i,n)
        l=l+1; % index for 1~31(or 30 or 29)
        k=k+1; % index for 1~366
        mth_d_txt(k,:)=[num2str(i+1988) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
        mth_d_txt_c{k,:}=[num2str(i+1988) '.' num2str(n,'%02d') '.'  num2str(j,'%02d')];
    end
    end
end

% pick matched date from water temp date
for i = 1:length(mth_d_txt_c)
       indx{i} = find([strcmp(mth_d_txt_c{i}, time_c)] == 1);
end
%%



%-- for var_name{ii}  처리부분
for ii = 1: length(var_name)
    %-- ii ==> 변수 순서 in var_name
    % variable_c=[r_txt_ud(2:end,9)];
    if( flag_ud ), variable_c=[r_txt_ud(1:end-1, var_col{ii})];
    else           variable_c=[r_txt_ud(2:end, var_col{ii})];  end
    
    for i = 1:length(variable_c)
        if strcmp(variable_c{i,1},'') == 1
           variable(i) = NaN;
        elseif strcmp(variable_c{i,1},'') == 0
           variable(i) = str2num(char(variable_c{i,1}));
        end
    end

    figure; hold on;
    plot(time, variable, 'linestyle','none','marker','.' ); xlim([x_tick(1) inf]); ylim([0 inf]);
    title( [fig_name,var_name{ii}] ,'fontsize',13);
    xlabel('date','fontsize',13); ylabel( [ var_name{ii}, ' (', dim_name{ii}, ')' ] ,'fontsize',13);
    grid on;
    set(gca,'fontsize',13,'fontweight','bold');
    %regression
    %slope y = b1*x
    nonan=variable(~isnan(variable));
    idx = find(isnan(variable) == 0);
    nanaxisx=find(isnan(variable) == 1);
    b1 = idx'\nonan'; % x\y for getting slop
    yCalc1 = b1*idx;
    % Slope & Intercept y = b0 + b1*x
    X = [ones(length(idx'),1) time(idx')]; %b_0 b_1 
    b = X\nonan';
    yCalc2 = X*b;
    plot(time(idx),yCalc2,'--','color','r','linew',2)
%     set( gca, 'xtick',x_tick,'xticklabel',datestr(x_tick),'fontsize',13); xticklabel_rotate
    xticklabel_rotate(x_tick',45, x_tick_c)
    gtext(['y = ',num2str(b(2),2),'x + ',num2str(b(1),2)],'Color','k','FontSize',16)
    print('-dpng', [fig_name,var_name{ii}]);
end

