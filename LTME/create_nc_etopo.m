clc; clear all; close all;
load('etopo_merged.mat');
depth_merged(find(depth_merged(:,:)<-10000))=depth_merged(find(depth_merged(:,:)<-10000))/1000.;
filename='merged_etopo.nc';
nccreate(filename,'lon','Dimensions',{'lon',7200});
nccreate(filename,'lat','Dimensions',{'lat',6000});
nccreate(filename,'depth_etopo','Dimensions',{'lon',7200,'lat',6000});
nccreate(filename,'depth_merged','Dimensions',{'lon',7200,'lat',6000});
ncwrite(filename,'lon',xx);
ncwrite(filename,'lat',yy);
ncwrite(filename,'depth_etopo',depth');
ncwrite(filename,'depth_merged',depth_merged');