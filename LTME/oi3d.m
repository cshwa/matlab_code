function [zi,err]=oi3d(x,y,t,z,xi,yi,Lx,Ly,Lt,E,Ln)
%function [zi,err]=oi3d(x,y,t,z,xi,yi,Lx,Ly,Lt,E,Ln)
%
% OI3D : Objective Mapping of a Scalar Field in 3D.
% 
%============================================================================
%
%function [zi,err]=oi3d(x,y,t,z,xi,yi,Lx,Ly,Lt,{E,Ln})
%
% DESCRIPTION :
%    Objetive mapping using Zmap(X)=Zobs*(inv(A)*C) where C is data-true
%    covariance vector and A is data-data covariance matrix. The covariance
%    function has the form of either the Gauss Ftn (F(r)=F0*exp(-r.^2/L^2))
%    or the Poisson Ftn (F(r)=F0*(1+r/L+r^2/3L^2)*exp(-r/L)) with the 
%    given decorrelation length scale (=L). Error estimation as a fraction
%    of true variance is calculated by Err/F0=1-C*inv(A)*C'.
%
% INPUT : x,y,t   = locations and time of observation(=z) in the Cartesian coordinate
%          z    = observation 
%         xi,yi = locations of mapped output(=zi) in the Cartesian coordinate
%          L    = decorrelation length scale of X,Y,T for the covariance function
%          E    = random or subgrid-scale noise as a fraction of the true 
%                 variance, e.g. small scale unresolved field, instrumental 
%                 uncertainty, or local sampling error. 
%                 (Optional, Default is E=0)      
%          Ln   = decorrelation length scale for the covariance function of
%                 the subgrid-scale noise. IF Ln=0, the covariance function
%                 has the form of delta function which describes an local
%                 instrumental uncertainty.  (Optional, Default is Ln=0)      
%
% OUTPUT : zi  = estimated values at (xi,yi)
%          err = estimated error relative as a fraction of true variance
%                (unknown)   
%
% Programmed by Dr. Y.O.Kwon
%
% Modified by Dr. Y.H. Kim
%
% REFERENCE : 
%    Bretherton, F.P., R.E. Davis and C.B. Fandry, 1976.
%    A techique for objective analysis and design of oceanographic experiments
%    applied to MODE-73. Deep-Sea Res., V23, pp559-582.
%
%++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


%------------------------------------------------------------------------
% Check the input arguments
%------------------------------------------------------------------------

  if nargin < 9
     error('yhkim_oi3d.m: not enough number of input arguments')
  elseif nargin < 10
     E=0; Ln=0;
  elseif nargin < 11
     Ln=0;
  end

  if any(size(x)-size(y))|any(size(y)-size(z))
     error('yhkim_oi3d.m: sizes of x,y,z must be same')
  end

  if any(size(xi)-size(yi))
     error('yhkim_oi3d.m: sizes of xi,yi must be same')
  end


%------------------------------------------------------------------------
% Set the variables
%------------------------------------------------------------------------
  x=x(:)'; y=y(:)'; t=t(:)'; z=z(:)';
  nobs=length(x);
  [nmapy,nmapx]=size(xi);
  nmap=nmapx*nmapy;
  xi=xi(:)'; yi=yi(:)';

%------------------------------------------------------------------------
% Defining Covariance Function (F=f(Lx,Ly,dx,dy))

%   F=inline('(1.0-(dx/Lx).^2-(dy/Ly).^2-(dt/Lt).^2).*exp(-(dx/Lx).^2-(dy/Ly).^2-(dt/Lt).^2)','Lx','Ly','Lt','dx','dy','dt');
   F=inline('exp(-(dx/Lx).^2-(dy/Ly).^2-(dt/Lt).^2)','Lx','Ly','Lt','dx','dy','dt');
%   F=inline('(1.0-(dx./Lx).^2-(dy./Ly).^2).*exp(-(dx./Lx).^2-(dy./Ly).^2)','Lx','Ly','Lt','dx','dy','dt');

%------------------------------------------------------------------------
% Calculate data-data Covariance matrix
%------------------------------------------------------------------------

  dx=abs(x(ones(1,nobs),:)'-x(ones(1,nobs),:));
  dy=abs(y(ones(1,nobs),:)'-y(ones(1,nobs),:));
  dt=abs(t(ones(1,nobs),:)'-t(ones(1,nobs),:));
  A=F(Lx,Ly,Lt,dx,dy,dt);
%  index = find(A<0);
%  A(index) = 0.0;
  A=A+E*eye(nobs);


%------------------------------------------------------------------------
% Calculate data-true Covariance matrix
%------------------------------------------------------------------------

  dx=abs(x(ones(1,nmap),:)'-xi(ones(1,nobs),:));
  dy=abs(y(ones(1,nmap),:)'-yi(ones(1,nobs),:));
  dt=abs(t(ones(1,nmap),:)');
  C=F(Lx,Ly,Lt,dx,dy,dt);
%  index = find(C<0);
%  C(index) = 0.0;  

%------------------------------------------------------------------------
% Calculate the estimated Field
%------------------------------------------------------------------------

  zi=z*(A\C);
  zi=reshape(zi,nmapx,nmapy);

%------------------------------------------------------------------------
% Calculate the error as a fraction of the true variance (unknown)
%------------------------------------------------------------------------

  err=diag(1-C'*(A\C));
  err=reshape(err,nmapx,nmapy);

return;
