function grd2 = gg(location)
  
% if nargin == 0
%   location = 'eas'; % default
% end
    scoord = [5 0.4 50 20]; 
    
switch location
  case { 'nena6','nena'}
    grd_file = '/home/wilkin/roms/nena/in/roms_nena_grid_6.nc';
    scoord = [5 0.4 50 30];
  case 'testk'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\kunsan\roms_grd_kunsan.nc' ;
    scoord = [5 0.4 50 20];
  case 'test'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_yw12_test2\roms_grd_yw12_ed.nc' ;
    scoord = [5 0.4 50 20];
  case 'eas'
    grd_file = 'd:\matlab\roms\input\roms_grd_jeju_sin2.nc'; %grd file path
    scoord = [5 0.4 50 20];
  case 'eas4'
    grd_file = 'D:\Roms\06-12-25(1_4)\roms_grid_4degree.nc'; %grd file path
    scoord = [5 0.4 50 20];
  case 'eas4sin'
    grd_file = 'd:\matlab\roms\input\roms_grd_jeju2.nc'; %grd file path
    scoord = [5 0.4 50 20];
 
    case 'yw_n10'
    grd_file = 'd:\matlab\roms\input\roms_grd_yw.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_n10_4m'
    grd_file = 'd:\matlab\roms\input\roms_grd_yw_4m.nc'; % new 1/10
    scoord = [5 0.4 50 20];
 
    
    
    case 'eas4_ghseo'
    grd_file = 'd:\matlab\roms\input\roms_grid_4degree_3.nc'; % new 1/10
    case 'eas4_ghseo1'
    grd_file = 'd:\matlab\roms\input\roms_grid_4degree_6.nc'; % new 1/10
    scoord = [5 0.4 50 20];    
    
     case 'eas12_ghseo'
%     grd_file = 'G:\ROMS\EAS10\roms_grid_0728_new.nc'; % new 1/10
    grd_file = 'G:\ROMS\EAS10\grd\roms_grid_2005.nc'; % new 1/10   
    scoord = [5 0.4 50 20]; 
    
    case 'yw_skku_curve_0530'
    grd_file = 'd:\matlab\roms\input\roms_grd_0530_curve_comb.nc'; % new 1/10
    scoord = [5 0.4 50 20];    
    case 'yw_skku_curve_1112'
    grd_file = 'd:\matlab\roms\input\roms_grd_1112_curve_comb.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_skku_curve_1113'
    grd_file = 'd:\matlab\roms\input\roms_grd_1113_curve_comb.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_skku_curve_1114'
    grd_file = 'd:\matlab\roms\input\roms_grd_1114_curve_comb.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_skku_curve_1121'
    grd_file = 'd:\matlab\roms\input\roms_grd_ywc1121.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_skku_curve_1122'
    grd_file = 'd:\matlab\roms\input\roms_grd_ywc1122.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    
    case 'yw_skku_curve_1001'
    grd_file = 'd:\matlab\roms\input\roms_grd_ywc1001.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    
    case 'yw_skkuc_2011'
%     grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_ys12\roms_grd_yw12comb.nc'; % new 1/10
    grd_file = 'd:\matlab\roms\input\roms_grd_yw12comb.nc'; % new 1/10
    scoord = [5 0.4 30 20]; 
    case 'yw_skkuc_2011r'
%     grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_ys12\roms_grd_yw12comb.nc'; % new 1/10
    grd_file = 'd:\matlab\roms\input\roms_grd_yw12comb3.nc'; % new 1/10
    scoord = [5 0.4 30 20]; 
    case 'yw_skkuc_2012'
%     grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_ys12\roms_grd_yw12comb.nc'; % new 1/10
    grd_file = 'd:\matlab\roms\input\roms_grd_yw12comb4.nc'; % new 1/10
    scoord = [5 0.4 30 20];  
    case 'ywr35_10'
    grd_file = 'E:\2011_ys12c_rot35_10\input\roms_grd_yw12comb4.nc'; % new 1/10
    scoord = [5 0.4 30 20];
    case 'yw12r35'    
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_yw12_rot35_4\roms_grd_yw12comb4.nc'; % new 1/10
    scoord = [5 0.4 30 20];
 
    
    
    case 'yw_e5_curve_0817'
    grd_file = 'd:\matlab\roms\input\roms_grd_0817_e5curve_comb.nc'; % new 1/10
    scoord = [5 0.4 50 20];  
     case 'yw_e5_curve_0818'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_etopo5_n10_curve\roms_grd_0818_e5curve_comb.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    
    
    




    
    
    
    case 'yw_test'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\roms_grd_yw_test.nc';
    scoord = [5 0.4 50 20];      
    
    case 'yw12Rc'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_yw12_rot\roms_grd_yw12c.nc';
    scoord = [5 0.4 50 20];  
    
    case 'yw12Rc18'
grd_file = 'G:\yw12c_r18\input\roms_grd_yw12comb_5m.nc';
    scoord = [5 0.4 50 20]; 
    
%     case 'yw12Rc25'
% grd_file = 'G:\yw12c_r25\input\roms_grd_yw12comb.nc';
%     scoord = [5 0.4 50 20]; 
    case 'yw12Rc25'
%     grd_file = 'F:\yw12Rc25\grd\roms_grd_12_r25.nc';
        grd_file = 'F:\yw12Rc25\grd\roms_grd_yw12comb.nc';
    scoord = [5 0.4 50 20];   
    
    case 'yw12Rc22'
%     grd_file =
%     'D:\matlab\roms_tool\make_grid\Build_grid\sin_yw12_rot22\roms_grd_yw12RC22.nc';
% grd_file = 'G:\yw12Rc22\input\roms_grd_yw12RC22_5m.nc';
grd_file = 'G:\yw12Rc22_2\grd\roms_grd_yw12comb.nc';

    case 'yw12Rc222'
grd_file = 'G:\yw12Rc22_2\grd2\roms_grd_yw12comb.nc';

    case 'yw12Rc12'
grd_file = 'd:\matlab\roms_tool\make_grid\Build_grid\sin_yw12_rot12\roms_grd_yw12comb_1.nc';
    case 'yw12Rc_r18'
grd_file = 'G:\yw12c_r18\input\roms_grd_yw12comb.nc';


%     case 'yw12Rc27new'
% grd_file = 'G:\yw12Rc27_3\input\roms_grd_12_r27_newE.nc';
%     case 'yw12Rc27newc'
% grd_file = 'G:\yw12Rc27_3\input\grd\roms_grd_yw12comb.nc';
%     case 'yw12Rc27test'
% grd_file = 'G:\yw12Rc27_3\input\grd\test.nc';
%     case 'yw12Rc27t1'
% grd_file = 'G:\yw12Rc27_3\input\grd\roms_grd_t1c.nc';
% 
    case 'yw12c'
grd_file = 'E:\2011_ys12c\input\roms_grd_yw12comb.nc';
grd_file = 'E:\2011_ys12c\input\roms_grd_yw12comb1.nc'; 
    scoord = [5 0.4 30 20]; 
    case 'yw12'
grd_file = 'E:\2011_ys12sms\input\roms_grd_yw12comb1.nc';
    scoord = [5 0.4 30 20]; 
    case 'yw12s'
grd_file = 'E:\2011_ys122s\input\roms_grd_yw12comb.nc';
grd_file = 'E:\2011_ys12_s2\input\roms_grd_yw12comb1.nc';
    scoord = [5 0.4 30 20]; 
    case 'yw12s1'
grd_file = 'E:\2011_ys12_s2\input1\roms_grd_yw12comb2.nc';
    case 'yw12s2'
grd_file = 'E:\2011_ys12_s2\input2\roms_grd_yw12comb2.nc';
    case 'yw12s3'
    scoord = [5 0.4 30 20]; 
grd_file = 'E:\2011_ys12_s2\input3\roms_grd_yw12comb2.nc';
    case 'yw12s4'
    scoord = [5 0.4 30 20]; 
grd_file = 'E:\2011_ys12_s2\input4\roms_grd_yw12comb2.nc';

    case 'yw122s'
grd_file = 'E:\2011_ys122s\input\roms_grd_yw12comb1.nc';
    scoord = [5 0.4 30 20]; 
    case 'yw122s2'
grd_file = 'E:\2011_ys122s\input\roms_grd_yw12comb2.nc';
    scoord = [5 0.4 30 20];     
    case 'yw12arc'
grd_file = 'E:\2011_ys12c_ysarc\grd\roms_grd_12_rot27.nc';
grd_file = 'E:\2011_ys12c_ysarc\grd\roms_grd_yw12comb2.nc';
    scoord = [5 0.4 30 20]; 
    case 'yw12arc1'
grd_file = 'E:\2011_ys12c_ysarc\grd\roms_grd_yw12comb21.nc';
    scoord = [5 0.4 30 20]; 
    
    case 'yw12Rc27'
grd_file = 'D:\matlab\roms\input\roms_grd_yw12comb27.nc';


    case 'yw12Rc27raw'
grd_file ='E:\2011_ys12c_rot27\grd_rot27re2\roms_grd_0530_curve_comb_ed_yzdep2.nc';
    scoord = [5 0.4 30 20]; 
    case 'yw12Rc27s'
grd_file ='E:\2011_ys12c_rot27\grd_rot27sms\roms_grd_yw12comb4.nc'; % new 1/10
    scoord = [5 0.4 30 20]; 
    

    case 'yw12Rc27_2'
grd_file = 'G:\yw12Rc27_2\input\roms_grd_yw12comb2.nc';
    case 'yw12Rc27_32'
grd_file ='G:\yw12Rc27_3\input\grd2\roms_grd_12_r27_new.nc'
grd_file ='G:\yw12Rc27_3\input\grd2\roms_grd_yw36comb.nc'


    case 'yw12Rc27_4'
grd_file = 'G:\yw12Rc27_2\input\roms_grd_yw12comb4.nc';
 
    case 'yw12Rc27_E5m2'
grd_file = 'G:\yw12Rc27_2\input\roms_grd_yw12combE_5m2.nc';
    case 'yw12Rc27_E5m2sk'
% grd_file = 'F:\04_ROMS_model_test\yw12Rc27_2\input\roms_grd_yw12combE_5m2sk.nc';
grd_file = 'E:\2011_ys12c_rot27\input\roms_grd_yw12combE_5m2sk.nc';
    scoord = [5 0.4 30 20];   
    
    case 'yw12Rc27_E5m2skE'
grd_file = 'G:\yw12Rc27_2\input\roms_grd_yw12combE_5m2skE.nc';

    case 'yw12Rc27_E5m3'
grd_file = 'G:\yw12Rc27_2\input\roms_grd_yw12combE_5m3.nc'; 
%     case 'yw12Rc27_new'
% grd_file = 'G:\yw12Rc27_2\input\roms_grd_yw12comb_sk.nc';  
 


    case 'yw12Rc30'
%     grd_file = 'F:\yw12Rc30\grd\roms_grd_12_r30.nc';
        grd_file = 'F:\yw12Rc30\grd\roms_grd_yw12comb.nc';
        grd_file = 'F:\yw12Rc30\grd\roms_grd_yw12combE3.nc';
    scoord = [5 0.4 50 20];   
    
    case 'yw3km'
%     grd_file = 'F:\yw12Rc30\grd\roms_grd_12_r30.nc';
        grd_file = 'H:\yw3km_my\input\roms_grd_ylw_expan5.nc';
    scoord = [5 0.4 50 20]; 
    
    
    case 'yw12Rc35'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_yw12_rot2_35\roms_grd_yw12Rc_35.nc';
    scoord = [5 0.4 50 20];    
    
    case 'yw12Rc1'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_yw12_rot\roms_eas_grid8.nc';
    scoord = [5 0.4 50 20]; 
    
    
    case 'yw12Rc_orth'
    grd_file='D:\matlab\roms_tool\etc_tools\Preprocessing_tools\roms_grd_12orth.nc';
    scoord = [5 0.4 50 20]; 
        case 'yw_test_tide'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin_test_tide\roms_grd_test_tide5.nc';
    scoord = [5 0.4 50 20];          
        case 'kodc20'
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\sin20\roms_grd_kodc20.nc'; % new 1/20
    scoord = [5 0.4 50 20];       
        case 'nwp'  %ccw
    grd_file = 'E:\test\ccw\roms_grd_1_4_smooth_hmax2.nc'; % new 1/20
    scoord = [5 0.4 50 20];   
    
        case 'nwps4'  % nwps4sin
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\NWPS4\roms_eas_grid8.nc';
    scoord = [5 0.4 50 20];   
    
        case 'nwp10new'  % nwps4sin
     grd_file = 'E:\2012_nwp10\input\roms_grid_10km_new.nc';      
    scoord = [5 0.4 30 20];  
        case 'nwp10new1'  % nwps4sin
     grd_file = 'E:\2012_nwp10\input\roms_grd_nwp10_1.nc';      
    scoord = [5 0.4 30 20];  
        case 'nwp10new2'  % nwps4sin
     grd_file = 'E:\2012_nwp2\grd\roms_grd_nwp10_2.nc';      
    scoord = [5 0.4 30 20];  
        case 'nwp10new3'  % nwps4sin
     grd_file = 'E:\2012_nwp2\grd2\roms_grd_nwp10_g2.nc';      
    scoord = [5 0.4 30 20];  
        case 'nwp10r'  % nwps4sin
     grd_file = 'E:\2012_nwp2\grd\roms_grd_nwp10r.nc';  
    scoord = [5 0.4 30 20];  
    
        case 'yw36'  % nwps4sin
%     grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\yw36_r35\roms_grd_yw36sin.nc';
%     grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\yw36_r35\roms_grd_yw36comb.nc';
    grd_file = 'G:\yw36Rc18_2\grd\roms_grd_yw36comb.nc';
    scoord = [5 0.4 50 20];  
    
        case 'yw36_3'  % nwps4sin
    grd_file =  'G:\yw36Rc18_3\grd\roms_grd_yw36comb.nc';
        case 'yw36_3_'  % nwps4sin
    grd_file =  'G:\yw36Rc18_3\grd\roms_grd_yw36comb1.nc';
    
        case 'yw36r12'  % nwps4sin
%     grd_file = 'G:\yw36Rc12\grd\roms_grd_yw36comb.nc';
    grd_file = 'G:\TEST\yw36Rc12\grd\use_sms\roms_grd_yw36combSW.nc'; 
    
        case 'yw36r12r'  % nwps4sin
%     grd_file = 'G:\yw36Rc12\grd\roms_grd_yw36comb.nc';
    grd_file = 'G:\yw36Rc12\grd\use_sms\roms_grd_yw36comb.nc';   
    
        case 'yw36r13'  % nwps4sin
    grd_file = 'G:\yw36Rc13\grd\roms_grd_yw36comb.nc';
    
        case 'yw36test'  % yw36test
%     grd_file = 'D:\matlab\roms_tool\etc_tools\Preprocessing_tools\roms_grd_yw36.nc';
    grd_file = 'D:\matlab\roms_tool\make_grid\Build_grid\yw36_r24\roms_grd_yw36comb.nc';
    grd_file = 'G:\yw36Rc18\grd\roms_grd_yw36comb.nc';
    scoord = [5 0.4 50 20];  
        case 'yw36orth'  % yw36test
    grd_file = 'G:\TEST\yw36orth\input\roms_grd_yw36comb_orth.nc';
    scoord = [5 0.4 50 20];  
        case 'yw36orth2'  % yw36test
    grd_file = 'G:\yw36orth2\input\roms_grd_yw36comb_orth2.nc';
    scoord = [5 0.4 50 20];
    
        case 'yw36kmt'  % yw36test
%     grd_file =     'H:\TEST\yw3km_orth3\grd\roms_grd_orth3_xyz.nc';
    grd_file =     'H:\yw3km_my\input\roms_grd_ylw_expan5.nc';
    
        case 'yw36orth3'  % yw36test
%     grd_file =     'H:\TEST\yw3km_orth3\grd\roms_grd_orth3_xyz.nc';
%     grd_file =     'H:\TEST\yw3km_orth3\grd\roms_grd_yw36orth3c.nc';   % bndy yw10c
 grd_file =     'H:\TEST\yw3km_orth3\grd\roms_grd_yw36orth3ec.nc';   % bndy yw10c
        case 'yw36orth32'  % yw36  with EAS10
    grd_file =     'H:\TEST\yw3km_orth3\grd\roms_grd_yw36orth3ec2.nc';  % bndy eas10
    scoord = [5 0.4 30 20]; 
    
        case 'yw36orth33'  % yw36  with EAS10
    grd_file =     'H:\TEST\yw3km_orth3\grd\roms_grd_yw36orth3ec3.nc';  % bndy eas10
        case 'yw36sq'  % yw36test
    grd_file = 'G:\yw36sq\grd\roms_grd_yw36comb_sq2.nc';
    scoord = [5 0.4 50 20];  
        case 'yw36ss'  % yw36test
%     grd_file = 'F:\yw36_ss\grd\roms_grd_yw36comb.nc';
%     grd_file = 'H:\yw3km_ss\grd\roms_grd_yw36comb2.nc';
    grd_file = 'H:\TEST\yw3km_ss\grd\roms_grd_yw36comb2.nc';
    scoord = [5 0.4 50 20];  
        case 'yw36ssSW'  % yw36test
%     grd_file = 'F:\yw36_ss\grd\roms_grd_yw36comb.nc';
    grd_file = 'F:\yw3km_ss\grd\roms_grd_yw36combSW.nc';
            case 'yw36ssSW2'  % yw36test
%     grd_file = 'F:\yw36_ss\grd\roms_grd_yw36comb.nc';
    grd_file = 'F:\yw3km_ss\grd\roms_grd_yw36combSW2.nc';
    
    scoord = [5 0.4 50 20];  
            case 'yw12ss'  % yw36test
    grd_file = 'F:\yw12_ss\grd\roms_grd_yw12comb.nc';
    scoord = [5 0.4 50 20];  
    case 'yw36r27'  % yw36test
    grd_file = 'F:\yw3km_27\grd\roms_grd_yw36comb.nc';
    scoord = [5 0.4 50 20];  
    case 'yw36r272'  % yw36test
    grd_file = 'H:\yw3km_27\grd\roms_grd_yw36comb2.nc';
    scoord = [5 0.4 50 20]; 
    case 'yw36r273'  % yw36test
    grd_file = 'H:\yw3km_27\grd\roms_grd_yw36comb3.nc';
    scoord = [5 0.4 50 20]; 
    case 'yw36r27nap'  % yw36test
    grd_file =  'H:\yw3km_27\grd\heavy_smoo\roms_grd_yw36c_sm4nap.nc';
%%
    case 'yw_skku_curve_0530edyz_dep'
    grd_file = 'd:\matlab\roms\input\roms_grd_0530_curve_comb_ed_yzdep.nc'; % new 1/10
    scoord = [5 0.4 50 20];  %%%    
    case 'yw_skku_curve_0530edyz'
    grd_file = 'd:\matlab\roms\input\roms_grd_0530_curve_comb_ed_yz.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_e5_curve_0817edyz'
    grd_file = 'd:\matlab\roms\input\roms_grd_0817_e5curve_comb_ed_yz.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_skc_0530edyz2'
    grd_file = 'd:\matlab\roms\input\roms_grd_0530_curve_comb_ed_yzdep2.nc'; % new 1/10
    scoord = [5 0.4 50 20]; 
    case 'yw_e5c_0817edyz2'
    grd_file = 'd:\matlab\roms\input\roms_grd_e5_curve_comb_ed_yzdep2.nc'; % new 1/10
    scoord = [5 0.4 50 20];
    
    case 'yw_skkuc1'
    grd_file = 'G:\ROMS\yw10c_0530_skku\roms_grd_0530_curve_comb1.nc'; % new 1/10
    scoord = [5 0.4 50 20];
 
    
    case 'yw_e5c_edyz3'
    grd_file = 'd:\matlab\roms\input\roms_grd_e5c_comb_edyzdep3.nc'; % new 1/10
    scoord = [5 0.4 50 20];    
    case 'yw_e5c_edyz3_30'
    grd_file = 'd:\matlab\roms\input\roms_grd_e5c_comb_edyzdep3.nc'; % new 1/10
    scoord = [5 0.4 50 30];    
    case 'yw_e5c_edyz3_40'
    grd_file = 'd:\matlab\roms\input\roms_grd_e5c_comb_edyzdep3.nc'; % new 1/10
    scoord = [5 0.4 50 40];  
    case 'ADD6'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_06_smooth_5.nc'; % new 1/10
    scoord = [3.5 0.2 10 40]; 
    case 'ADD4'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_04_smooth_2.nc'; % new 1/10
    scoord = [3.5 0.2 10 40]; 
    case 'ADD'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_03.nc'; % new 1/10
    scoord = [3.5 0.2 10 40]; 
    case 'ADD5'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_05_smooth_3.nc'; % new 1/10
    scoord = [3.5 0.2 10 40]; 
    case 'ADD7'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_07_smooth_5.nc'; % new 1/10
    scoord = [3.5 0.2 10 40]; 

    case 'ADD6_ori'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_06.nc'; % new 1/10
    scoord = [3.5 0.2 10 40]; 

    case 'ADD8'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_08.nc'; % new 1/10
    scoord = [3.5 0.2 10 40];     
    
    case 'ADD10'
    grd_file = 'D:\add_ini_bry_grd\grid\roms_grid_ADD_10.nc'; % new 1/10
    scoord = [3.5 0.2 10 40];        
    
    case 'ADD2_10'
    grd_file = 'D:\add2_ini_bry_grd\grid\roms_grid2_ADD_10.nc'; % new 1/10
    scoord = [3.5 0.2 10 40];
       case 'ADD2_10ep'
    grd_file = 'D:\add2_ini_bry_grd\grid\roms_grid2_ADD_10_ep.nc'; % new 1/10
    scoord = [3.5 0.2 10 40];  
           case 'ADD2_08ep'
    grd_file = 'D:\add2_ini_bry_grd\grid\roms_grid2_ADD_08_ep.nc'; % new 1/10
    scoord = [3.5 0.2 10 40];  
           case 'ADD2_04_2epo'
    grd_file = 'D:\add2_ini_bry_grd\grid\roms_grid2_ADD_04_2_epo.nc'; % new 1/10
    scoord = [5 0.4 5 40];     
           case 'ADD2_08_2ep'
    grd_file = 'D:\add2_ini_bry_grd\grid\roms_grid2_ADD_08_2_ep.nc'; % new 1/10
    scoord = [5 0.4 5 40];  
            case 'ADD3_04e'
    grd_file = 'D:\add3_ini_bry_grd\grid\roms_grid3_ADD_04_e.nc'; % new 1/10
    scoord = [3.5 0.2 10 40];  
            case 'ADD3_04_2e'
    grd_file = 'D:\add3_ini_bry_grd\grid\roms_grid3_ADD_04_2_e.nc'; % new 1/10
    scoord = [5 0.4 5 40];
                case 'ADD3_07_2e'
    grd_file = 'D:\add3_ini_bry_grd\grid\roms_grid3_ADD_07_2_e.nc'; % new 1/10
    scoord = [5 0.4 5 40];  
           case 'ADD3_08_2ep'
    grd_file = 'D:\add3_ini_bry_grd\grid\roms_grid3_ADD_08_2_ep.nc'; % new 1/10
    scoord = [5 0.4 5 40];  
           case 'ADD4_05_2ep'
    grd_file = 'D:\add4_ini_bry_grd\grid\roms_grid4_ADD_05_2_ep.nc'; % new 1/10
    scoord = [5 0.4 5 40];  
                 case 'ADD3_10_2e'
    grd_file = 'D:\add3_ini_bry_grd\grid\roms_grid3_ADD_10_2_e.nc'; % new 1/10
    scoord = [5 0.4 5 40];  
                case 'reanal_25km'
    grd_file = 'D:\monthly_25km\roms_grid_4degree_3.nc'; % new 1/10
    scoord = [5 0.4 5 20];  
                    case 'reanal_10km'
    grd_file = 'D:\monthly_10km\roms_grid02_2.nc'; % new 1/10
    scoord = [5 0.4 5 20];  
    
    
end


disp(' ')
disp([ 'Loading ROMS grd for application: ' location])
disp([ 'using grid file ' grd_file])
disp(' ')
grd2 = roms_get_grid2(grd_file,scoord);

  
