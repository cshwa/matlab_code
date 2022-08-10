function [var_rho]=w2rho_2d(var_w)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) 2000 IRD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
% transfert a field at u points to a field at rho points
%
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[Mp,L]=size(var_w);
Mm=Mp-1;
var_rho=zeros(Mm,L);
var_rho(1:Mm,:)=0.5*(var_w(1:Mm,:)+var_w(2:Mp,:));
% var_rho(:,1)=var_rho(:,2);
% var_rho(:,Lp)=var_rho(:,L);
return

