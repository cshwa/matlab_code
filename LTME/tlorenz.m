function [xp,yp,zp] = tlorenz(x,y,z,xb,yb,zb,p,r,b)

% function [xp,yp,zp] = tlorenz(x,y,z,xb,yb,zb,p,r,b)
%
% Tangent Linear Lorenz Model
% Programmed by Y.H. Kim for the summer school on Aug. 2007 
%
% x, y, z    : input k-th results of state vector
% xb, yb, zb    : input k-th results of nonlinear lorenz model
% p, r, b    : input parameter
% xp, yp, zp : k+1-th results of state vector
%

xp = -p*x + p*y;
yp = (r-zb)*x - y - xb*z;
zp = yb*x + xb*y - b*z;
end