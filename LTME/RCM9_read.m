cc
% load './data_fin/EC1_SBE_combine.mat'
fpath='D:\work\exploration\2012_Jul_Eardo_EC1_recover_mooring\data\'

fname='RCM9_1241_DSU_15194_EC1_1000m_2011Mar_2012Jul.asc';

% [2011,3,12,16,30,0;] [2012,7,28,8,0,0;]
rd=importdata([fpath fname]);
rdd=rd;
idd=find(rdd(:,8)==400);
dataa=rdd(idd,:);

data=dataa(12:24187,:);
dataa=data;
time=[dataa(:,4) dataa(:,3) dataa(:,2) dataa(:,5) dataa(:,6) dataa(:,7) ];
jjd=datenum([2011 3 12 16 30 0]):1/24/2:datenum([2012 7 28 8 0 0]);


time=datevec(jjd);
jd=datenum(time);


 itmp=data(:,11);
icond=data(:,12);
ipre=data(:,13);
 idr=data(:,10);
 isp=data(:,9);
 
 %temp
 	ta = -3.098E00;
	tb =  8.960E-3;
	tc = -3.476E-7;
	td =  1.134E-10 ;
    temp = ta + (tb*itmp) + (tc*itmp.^2) + (td*itmp.^3);
    


%direction
	da = 0                                                                       %direction
	                                    
	db = 3.516E-1  
    dir0 = da + db*idr                                   %direction
	dir = dir0
%spd
sa = 0
	sb = 2.933E-1                                              %speed 단위가 cm/s
	spd = sa + sb*isp
%pre
    pa = -5.534E-1;
    pb = 2.043E-2;
    pc = 2.262E-7;
    pd = 0
      pre0 = pa + (pb*ipre) +(pc*ipre.^2)                                      % Mpa = 10^6pa
	pre = pre0  * 100  % Mpa  -> dbar
    
    %con -> sal
%     SBE37 Sal 계산

% S37Prs01 = sw_pres(S37Dpt01,GPSlat);
% sal = sw_salt(cond./sw_c3515(), temp,pre);
% [ur vr]=wind2cur(spd,dir,0);
% u=-ur;v=-vr;
% figure
% plot(u); hold on; plot(v,'r')
% [u_mag v_mag]=magvari(u,v,-8)

save 'RCM9_1241_DSU_15194_EC1_1000m_2011Mar_2012Jul.mat' pre spd dir  u_mag v_mag temp time