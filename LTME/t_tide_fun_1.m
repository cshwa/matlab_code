function [pout]=t_tide_fun(raw_data,inter,lat,stat_t,output,index_sv)
% Setup for t_tide run.
%  
% raw_data
% inter = interval
% lat = latitude (decimal degrees (+north))
% stat_t = start time (decimal day (matlab DATENUM scalar))
% output = result output
% index_sv = choose scaler data(1) or vector data(2)
%
% By Kim Dae Hyun
% Ver 0.98 Test run.
% Ver 0.99 made Function 2002. 1. 23


if(index_sv==0)
	% Define inference parameters for Scalar data.
     infername=['P1';'K2'];
     inferfrom=['K1';'S2'];
     infamp=[0.33093;0.27215];
     infphase=[-7.07;-22.40];
elseif(index_sv==1)
	% Define inference parameters for Vector data.
     infername=['P1';'K2'];
     inferfrom=['K1';'S2'];
     infamp=[0.311807 0.197553;0.191983 0.212745];
     infphase=[-4.5 -0.2;-14.8 27.2];         
end

[nameu,fu,tidestruc,pout]=t_tide(raw_data,...
	'interval',inter,...
	'inference',infername,inferfrom,infamp,infphase,...
	'latitude',lat,...              
	'start',stat_t,...	
	'secular','mean',...
	'output',output,...
     'shallow','M10',...                
     'error','cboot',...                   
     'synthesis',0);

% Z0 include
%	pout = pout + z0;        

