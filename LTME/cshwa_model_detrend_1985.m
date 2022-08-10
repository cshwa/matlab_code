% T_DEMO - demonstration of capabilities.
% Short example of capabilities of tidal analysis toolbox.
%
% In this example, we 
%         a) do nodal corrections for satellites, 
%         b) use inference for P1 and K2, and
%         c) force a fit to a shallow-water constituent.

% Version 1.0

% Sept/2014 - it was pointed out that the elevations are in feet, not meters.
%             plot labels changed. Also added a note about time zones.

clear all; clc; close all; 
echo on
       echo on
       % Load the example.
         load cshwa_1985_time_fix.mat  % this is a local time
%          time(4545:end)=[];
         load('1985_regrid_v3_zeta.mat');
         total_temp(length(time)+1:end)=[]; %model has more
         ssh=total_temp;
         tg_h=(ssh-mean(ssh)).*100;
         elevation = tg_h .*0.0328084;  %cm to feet
         tg_h_f = tg_h .*0.01; % cm to m
         %%% determin when do start the detrend
         tt=datestr(time,'mm/dd'); tt1=datestr(time,'mm/dd HH'); %% 03/06 22
         plot(elevation);hold on; % line(1602,0:0.001:max(elevation),'color','r');
         elevation1 = elevation;
%          elevation1 = nandetrend(elevation);
         % case GMT, convert to GMT (local = GMT - 9.00) so we add 9 hour
%          plot(nandetrend(elevation1));
         time_f=time;
         tt_f=datestr(time_f,'mm/dd HH'); %% tt_f(356,:) >>> 03/23 21, tt_f(1965,:) >>> 05/29 22
         
         %%% plot trend %%
         x=(1:length(tg_h_f));
         idx = find(isnan(tg_h_f)==1)
         p1 = polyfit(x, tg_h_f,1);
         y1 = polyval(p1,x');
         
        figure;
        plot(tg_h.*0.01,'linew',2);hold on;
        plot(1:length(tg_h),y1,'r','linew',2);
        title('TG - detrend');
        set(gca,'linew',2);
        xlabel('Time (date)');
        ylabel('elevation (m)');
        set(gca,'xtick',1:240:length(time));
        set(gca,'xticklabel',tt(1:240:length(time),:));
%         set(gca,'XTickLabelRotation',45);
        set(gca,'fontsize',20);
         
%%%%%%%%%%%did in the 2016 matlab%%%%%%%%%%%%%%%%%%%%%
%        t1 = datetime(2018,03,23,12,00,00); >> 21,00,00 (GMT)
%        t2 = datetime(2018,05,29,13,30,00); >> 22,00,00  (GMT)
%        t = t1:hours(1):t2
%        time = datenum(t)
%        save('cshwa_TD_time_fix','time')
%
%      cf.)  time = datenum('01/01/2015 00:00:00','mm/dd/yyyy HH:MM:SS')
%            daysact('01-01-1700','01-01-2014') %% conform days from 1700 is right. 
%            leapyear(2014) %% conforming year 365 or 366
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
       % Define inference parameters.
       infername=[];
       inferfrom=[];
       infamp=[];
       infphase=[];
%        infername=['P1';'K2'];
%        inferfrom=['K1';'S2'];
%        infamp=[.33093;.27215];
%        infphase=[-7.07;-22.40];
       
       % The call (see t_demo code for details).
       [tidestruc,pout]=t_tide(elevation1,...
       'interval',1, ...                     % hourly data
       'start',time_f(1),...               % start time is datestr(tuk_time(1))
       'latitude',34.74721840,...               % Latitude of obs
       'inference',infername,inferfrom,infamp,infphase,...     % Add a shallow-water constituent 
       'error','linear',...                   % coloured boostrap CI
       'synthesis',1);                       % Use SNR=1 for synthesis. 

     % Note - in the demo I use tuk_time which is MST. Strictly speaking,
     % in order to get Greenwich phase one should use Greenwich TIME as well,
     % which means adding 6 hours.
     % However, I stuck with doing things the sloppy not-quite-correct way so that 
     % results would comparable to those published in Mike Foreman's data report on the Fortran code.
     
       echo off

%        pout=t_predic(tuk_time,tidestruc,,...
%                      'latitude',69+27/60,...
%                      'anallength','full',...
%                      'synthesis',1);

% meter = ft/3.281

csh_name = {'M2', 'S2', 'N2', 'K2', 'K1', 'O1','P1','Q1'};
for i = 1:8
       srt_ind(i)=find(strcmp(csh_name(i),tidestruc.name)==1);
       n_srt_ind(:,i)=strcmp(csh_name(i),tidestruc.name);
end
no_G8_ind=find(sum(n_srt_ind,2)==0);

tidestruc.name(srt_ind,:) %confirm G8 contituents
tidestruc_8 = struct('name', tidestruc.name(srt_ind,:), ...
    'freq', tidestruc.freq(srt_ind,:), 'tidecon',tidestruc.tidecon(srt_ind,:), 'type',tidestruc.type);

tidestruc_no_G8 = struct('name', tidestruc.name(no_G8_ind,:), ...
    'freq', tidestruc.freq(no_G8_ind,:), 'tidecon',tidestruc.tidecon(no_G8_ind,:), 'type',tidestruc.type);

       pout_8=t_predic(time_f,tidestruc_8,...
                     'latitude',34.74721840,...
                     'anallength','nodal',...
                     'synthesis',0);
                 
       pout_no_G8_0=t_predic(time_f,tidestruc_no_G8,...
                     'latitude',34.74721840,...
                     'anallength','nodal',...
                     'synthesis',0);
                 
      pout_no_G8=t_predic(time_f,tidestruc_no_G8,...
                     'latitude',34.74721840,...
                     'anallength','nodal',...
                     'synthesis',1);
      
     pout_raw=t_predic(time_f,tidestruc,...
                     'latitude',34.74721840,...
                     'anallength','nodal',...
                     'synthesis',1);           
                 
                 
 save('MODEL_1985_G8_nondet_fix.mat','pout_*','tidestruc*');
 return
 tidestruc.name(~srt_ind,:) %confirm non-G8 contituents
 

% save('model_1980_v4.mat')

clf;orient tall;
subplot(411);
% plot(1:1610,[elevation1 pout]);
plot(1:length(time_f), pout./3.281);
% line(time-datenum(2018,1,0),elevation1-pout,'linewi',2,'color','r');
xlabel('Days in 1980','fontweight','bold');
ylabel('Elevation (m)','fontweight','bold');
% set(gca,'ytick',-0.8:0.2:0.8);
% set(gca,'yticklabel',-0.8:0.2:0.8);
set(gca,'xtick',1:600:length(time_f));
set(gca,'xticklabel',tt_f(1:600:length(time_f),:));
% text(190,5.5,'Original Time series','color','b');
% text(190,4.75,'Tidal prediction from Analysis','color',[0 .5 0]);
% text(190,4.0,'Original time series minus Prediction','color','r');
% legend({'Original data','Tidal prediction', 'diff'})
legend({'Tidal prediction'});grid('on');
title('TG-yeosu elevation');  xlim([1 4580]); % ylim([-0.8 0.8]);
set(gca,'fontweight','bold');

subplot(412);
fsig=tidestruc.tidecon(:,1)./3.281>tidestruc.tidecon(:,2)./3.281; % Significant peaks
semilogy([tidestruc.freq(~fsig),tidestruc.freq(~fsig)]',[.0005*ones(sum(~fsig),1)./3.281,tidestruc.tidecon(~fsig,1)./3.281]','.-r');
line([tidestruc.freq(fsig),tidestruc.freq(fsig)]',[.0005*ones(sum(fsig),1)./3.281,tidestruc.tidecon(fsig,1)./3.281]','marker','.','color','b');
line(tidestruc.freq,tidestruc.tidecon(:,2)./3.281,'linestyle',':','color',[0 .5 0]);
set(gca,'ylim',[.0005 1],'xlim',[0 .5]);
xlabel('frequency (cph)','fontweight','bold');
text(tidestruc.freq,tidestruc.tidecon(:,1)./3.281,tidestruc.name,'rotation',60,'vertical','base');
ylabel('Amplitude (m)','fontweight','bold');
text(.27,.4,'Analyzed lines with 95% significance level');
text(.25,.2,'Significant Constituents','color','b');
text(.25,.1,'Insignificant Constituents','color','r');
text(.25,.05,'95% Significance Level','color',[0 .5 0]);
xlim([0 .33]);
set(gca,'fontweight','bold');

subplot(413);
errorbar(tidestruc.freq(~fsig),tidestruc.tidecon(~fsig,3),tidestruc.tidecon(~fsig,4),'.r');
hold on;
errorbar(tidestruc.freq(fsig),tidestruc.tidecon(fsig,3),tidestruc.tidecon(fsig,4),'o');
hold off;
set(gca,'ylim',[-45 360+45],'xlim',[0 .33],'ytick',[0:90:360]);
xlabel('frequency (cph)','fontweight','bold');
ylabel('Greenwich Phase (deg)','fontweight','bold');
text(.27,330,'Analyzed Phase angles with 95% CI');
text(.25,290,'Significant Constituents','color','b');
text(.25,250,'Insignificant Constituents','color','r');
set(gca,'fontweight','bold');

subplot(414);
ysig=elevation1;
yerr=elevation1-pout;
nfft=389;
bd=isnan(ysig);
gd=find(~bd);
bd([1:(min(gd)-1) (max(gd)+1):end])=0;
ysig(bd)=interp1(gd,ysig(gd),find(bd)); 
%[Pxs,F]=psd(ysig(isfinite(ysig)),nfft,1,[],ceil(nfft/2));
[Pxs,F]=pwelch(ysig(isfinite(ysig)),hanning(nfft),ceil(nfft/2),nfft,1);
Pxs=Pxs/2;
%%[Pxso,Fo]=psd(ysig(isfinite(ysig)),nfft,1,[],ceil(nfft/2));

%[Pxs,F]=pmtm(ysig(isfinite(ysig)),4,4096,1);
yerr(bd)=interp1(gd,yerr(gd),find(bd)); 
%[Pxe,F]=psd(yerr(isfinite(ysig)),nfft,1,[],ceil(nfft/2));
[Pxe,F]=pwelch(yerr(isfinite(ysig)),hanning(nfft),ceil(nfft/2),nfft,1);
Pxe=Pxe/2;
%[Pxe,F]=pmtm(yerr(isfinite(ysig)),4,4096,1);

semilogy(F,Pxs);
line(F,Pxe,'color','r');
xlabel('frequency (cph)','fontweight','bold');
ylabel('m^2/cph','fontweight','bold');
text(.17,1e4,'Spectral Estimates before and after removal of tidal energy');
text(.25,1e3,'Original (interpolated) series','color','b');
text(.25,1e1,'Analyzed Non-tidal Energy','color','r');
xlim([0 .33]);
ylim([10^-5 10^10]);
set(gca,'fontweight','bold');



