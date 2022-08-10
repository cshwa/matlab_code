% T_DEMO - demonstration of capabilities.
% Short example of capabilities of tidal analysis toolbox.
%
% In this example, we 
%         a) do nodal corrections for satellites, 
%         b) use inference for P1 and K2, and
%         c) force a fit to a shallow-water constituent.

% Version 1.0
     
       echo on
       % Load the example.
       load jeju.dat %t_exajejule
       input=jeju(:,1);
       out=['jeju_u.dat']
      
      % The call (see t_demo code for details).
       t_eleu=input;
       [tidestruc,pout]=t_tide(t_eleu,...
       'interval',3/6, ...                   % hourly data
       'start',[2004,5,27,11,00,0],...        % start time is datestr(tuk_time(1))
       'latitude',33+12/60+4/3600,...       % Latitude of obs
       'output',out,...
       'shallow','M10',...                   % Add a shallow-water constituent 
       'error','cboot',...                   % coloured boostrap CI
       'synthesis',1);                       % Use SNR=1 for synthesis. 


       echo off
       
       
              echo on
       % Load the exajejule.
       load jeju.dat %t_exajejule
       input=jeju(:,2);
       out=['jeju_v.dat']
      
       % The call (see t_demo code for details).
       t_elev=input;
       [tidestruc,pout]=t_tide(t_elev,...
       'interval',3/6, ...                   % hourly data
       'start',[2004,3,26,15,00,0],...        % start time is datestr(tuk_time(1))
       'latitude',33+12/60+4/3600,...       % Latitude of obs       
       'output',out,...
       'shallow','M10',...                   % Add a shallow-water constituent 
       'error','cboot',...                   % coloured boostrap CI
       'synthesis',1);                       % Use SNR=1 for synthesis. 


       echo off


       echo on
       % Load the exajejule.
       load jeju.dat;  
       
       u2x=jeju(:,1);       
       v2y=jeju(:,2)*i;
       uv2xy=u2x+v2y;
       data(:,1)=uv2xy;
       out=['jeju_mm.dat']     
       % Define inference parameters.    
       % The call (see t_demo code for details).
       t_ele=data;
       [tidestruc,pout]=t_tide(t_ele,...
       'interval',3/6, ...                   % hourly data
       'start',[2004,3,26,15,00,0],...        % start time is datestr(tuk_time(1))
       'latitude',33+12/60+4/3600,...       % Latitude of obs       
       'output',out,...
       'shallow','M10',... 
       'error','cboot',...                   % coloured boostrap CI
       'synthesis',1);                       % Use SNR=1 for synthesis. 


       echo off

    %    pout=t_predic(tuk_time,tidestruc,,...
    %                  'latitude',69+27/60,...
    %                  'synthesis',1);

