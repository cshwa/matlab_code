function [xp,yp,zp] = florenz(x,y,z,p,r,b)

% function [xp,yp,zp] = florenz(x,y,z,p,r,b)
%
% Nonlinear Lorenz Model
% Programmed by Y.H. Kim for the summer school on Aug. 2007 
%
% x, y, z    : input k-th results of state vector
% p, r, b    : input parameter
% xp, yp, zp : k+1-th results of state vector
%

xp = -p*x + p*y;
yp = x*(r-z) - y;
zp = x*y - b*z;
end