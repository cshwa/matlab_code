close all; clear; clc; 
load uv.mat

%meanwd=atand(mean(sind(wd))/mean(cosd(wd)));
%meanwd=mod(meanwd,360)
wtime=[6/24:6/24:365];
yyc=0;
wm=zeros(386,488,10);wd=zeros(386,488,10);
f_u=zeros(386,488,1460);f_v=zeros(386,488,1460);
for tt=1:1:length(wtime)
    for yy=styear:endyear
        yyc=yyc+1;
        filenameu=['auto_ERA5_',num2str(yy),'_Uwind.nc'];
        filenamev=['auto_ERA5_',num2str(yy),'_Vwind.nc'];
        Uwind=ncread(filenameu,'Uwind',[1 1 tt],[Inf Inf 1]);
        Vwind=ncread(filenamev,'Vwind',[1 1 tt],[Inf Inf 1]);
        [twm,twd]=uv2compass(Uwind,Vwind);
        wm(:,:,yyc)=twm;wd(:,:,yyc)=twd;
    end
    f_tu=mean(wm,3).*mean(sin(wd*pi/180),3);
    f_tv=mean(wm,3).*mean(cos(wd*pi/180),3);
    f_u(:,:,tt)=f_tu;f_v(:,:,tt)=f_tv;yyc=0;
end

[twm,twd]=uv2compass(u,v);

f_tu=mean(twm,2).*mean(sin(twd*pi/180),2);
f_tv=mean(twm,2).*mean(cos(twd*pi/180),2);

smean=mean(sind(twd),2);
cmean=mean(cosd(twd),2);

phi_mean=atand(smean / cmean);

if cmean < 0
    phi_mean = phi_mean + 180;
elseif smean < 0 
    phi_mean = phi_mean + 360;
end

f_tu=mean(twm,2) * sind(phi_mean);
f_tv=mean(twm,2) * cosd(phi_mean); 

