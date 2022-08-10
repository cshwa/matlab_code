function [xap,yap,zap] = alorenz(xa,ya,za,xb,yb,zb,p,r,b)

% function [xap,yap,zap] = alorenz(xa,ya,za,xb,yb,zb,p,r,b)
%
% Adjoint Lorenz Model
% Programmed by Y.H. Kim for the summer school on Aug. 2007 
%
% xa, ya, za    : input k+1-th results of state vector
% xb, yb, zb    : input k-th results of nonlinear lorenz model
% p, r, b    : input parameter
% xap, yap, zap : k-th results of state vector
%

xap = -p*xa + (r-zb)*ya + yb*za;
yap = p*xa - ya + xb*za;
zap = -xb*ya - b*za;
end