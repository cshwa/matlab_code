% bfgs.m -- BFGS/Armijo/Skipping script for unconstrained minimization
%
% Script applies lightly-modified Newton method to a given function, 
% using workspace variables described below:
%   fname  - string containing name of function 
%   dfname - string containing name of gradient function
%   x0     - column vector containing starting point
% Stopping test parameters are defined near the top, and can be
% changed without harming what remains.  These are
%   eta1   - relative error among components of the input vector x
%   eta2   - tolerance for the gradient being nearly zero
%   Kmax   - maximum number of iterations

% (c) 2001, Philip D. Loewen
% 22 Feb 01  --  Original
% 23 Feb 01  --  Cosmetics
% 24 Feb 01  --  Reformatted final report

%
% Aug. 2007, Y.H. Kim
% For the summer school, solution of parameter estimation.
%

error_ratio = 1.30;

fname = 'func';
dfname = 'dfunc';

x0 = [1; 3; 5];  % Initial value of state vector
p0 = 10.0;       % Initial value of parameter p
r0 = 32.0;       % Initial value of parameter r
b0 = 2.66666667; % Initial value of parameter b

x0 = error_ratio * x0; % Add error to the state vector
pc = error_ratio * p0; % Add error to the parameter p
rc = r0;
bc = b0;

eta1 = 100*sqrt(eps)
eta2 = 10000*eps
Kmax = 100;       % Students were asked to use 25

Zeps = eps^(3/4); % Prevents division by 0 in stopping tests.

% Armijo line search parameters
alpha = 1.0e-3;  % Slope multiplier for "sufficient decrease" condition
beta  = 0.9;     % Geometric factor for backtracking iteration

% Set up table titles
disp([' ']);
disp(['BFGS/Armijo/Skipping minimization of function "',fname,'",']);
disp(['using gradient "',dfname,'":']);

xstring = '';
for jj=1:length(x0),
  xstring = [xstring,'x_k(',int2str(jj),')        '];
end;
disp([' ']);
disp(['  k     ',xstring,'f(x_k)        e_1           e_2         lambda']);

% Start Clock
tic;

% Set current point, function, and gradient. Define identity matrix.
xc = [x0; pc; rc; bc];
fc = eval([ fname,'(xc)']);
[gc] = eval([dfname,'(xc)']);
Id = eye(length(xc));

% Approximate Hessian Hc; Approximate Inverse Hessian Wc
Hc = abs(fc)*Id;
Wc = 1/abs(fc)*Id;

% Build output line
disp(['  0', sprintf('  %12.4e',[xc',fc])]);

index = 0;

for k=1:Kmax,

% Debugging line compares approx Hessian and approx inverse Hessian
% (Need to uncomment several lines below to evaluate/update Hc, sH.)
% disp([sprintf('%3d',k),sprintf('  norm(Hc*Wc-I)=%8.5f,',norm(Hc*Wc-Id,2))]);

%  sH = -Hc \ gc'; % QuasiNewton step computed using Hessian
  sW = -Wc * gc'; % QuasiNewton step computed using inverse-Hessian

  s  = sW;
  
  % Test if s points in a descent direction.
  if (gc*s > -eps^0.25*norm(gc,2)*norm(s,2)),
    % Newton step is rather poor descent direction.
    %  This shouldn't happen.
    disp([sprintf('%3d',k),'  trouble:  ',...
         's  is a poor descent direction.   Continuing ...']);
  end;

  % Take an Armijo step in direction s.
  m  = gc*s;       % Negative slope in direction s from point xc
  lam = 1;         % Stretch factor lambda
  while 1,         % Repeat forever:
    xn = xc + lam*s;                                 % New point
    fn = eval([ fname,'(xn)']);                      % Function value
    if ( (fn-fc) <= (alpha*lam*m)+eps ), break, end; % Test acceptability
    lam = lam*beta;                                  % Shrink lambda, retry
  end;
  % Get here with new point xn and fcn value fn in good shape.
  % Work out new gradient too.
  gn = eval([dfname,'(xn)']); 
  dx = xn - xc;   % Record change in x-vector.
  dg = gn - gc;   % Record change in gradient.
  % Evaluate stopping criteria.
  e1 = max(abs(dx') ./ max([abs(xc');Zeps*ones(size(xc'))])); % Rel chg in x
  e2 = max( (abs(gn') .* abs(xn)) / max([abs(fn),Zeps]) );  % Max rel chg in df

  % Build output line
  ol = sprintf('%3d',k);
  ol = [ol, sprintf('  %12.4e',xn')];
  ol = [ol, sprintf('  %12.4e',[fn,e1,e2])];
  ol = [ol, sprintf('  %8.6f',lam)];
  disp(ol);

  xc = xn;                     % Update the x-vector
  fc = fn;                     % Update the function value
  gc = gn;                     % Update the gradient

  index = index + 1;
  bfgs_ana(index).x = xc;
  bfgs_ana(index).cost = fc;
  
  % Act on stopping criteria.
  if (e1<eta1)|(e2<eta2), break, end
  
  % Approximate Hessian Update
  yts = dg*dx;
  if yts > sqrt(eps)*norm(dg,2)*norm(dx,2),
    y = dg'; s = dx;

%   Hcs = Hc*s;                                   %% Use with H
%   Hc = Hc + (y*y')/yts - (Hcs*Hcs')/(s'*Hcs);   %% Use with H

    Wy = Wc*y;
    Wyst = Wy*s';
    Wc = Wc + (yts+y'*Wy)/(yts)^2*s*s' - (Wyst' + Wyst)/yts;
  else
    disp([sprintf('%3d',k),'  ** Note:  ',...
         's''y is very small.  Hessian update skipped.']);
  end;

end;    % Finished main iteration.

tot_time = toc;

% Print results.
disp([' ']);
disp(['Best point:  x'' =',sprintf(' %22.16e',xc),'.']);
disp(['Best value:  ',fname,'(x) =',sprintf(' %22.16e',fc),'.']);
disp(['Step count:  ',int2str(k),' BFGS quasi-Newton steps.']);

ol = '';
if e1<eta1, 
  ol = 'e1 < eta1'; 
  if e2<eta2,
    ol = [ol,' and '];
  end;
end;
if e2<eta2,
  ol = [ol,'e2 < eta2'];
end; 
if k==Kmax,
   ol = 'iteration limit reached';
end;
disp(['Stopped by:  ',ol,'.']);

disp(['Elapsed time:  ',int2str(tot_time),' total--average ',int2str(tot_time/k),' per step.']);
disp(' ');

save bfgs_ana.mat bfgs_ana;