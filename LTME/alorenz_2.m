function [xap,yap,zap,pep,rep,bep] = alorenz(xa,ya,za,xb,yb,zb,p,r,b)

% function [xap,yap,zap,pep,rep,bep] = alorenz(xa,ya,za,xb,yb,zb,p,r,b)
%
% adjoint Lorenz Model
% Programmed by Y.H. Kim for the summer school on Aug. 2007 
%
% xa, ya, za : input k+1-th backward results of state vector
% xb, yb, zb : input background state vector
% p, r, b    : input k+1-th backward results of parameter
% xap, yap, zap : k-th backward results of state vector
% pep, rep, bep : k-th backward results of parameter
%

xap = -p*xa + (r-zb)*ya + yb*za;
yap = p*xa - ya + xb*za;
zap = -xb*ya - b*za;
pep = (yb-xb)*xa;
rep = xb*ya;
bep = -zb*za;
end