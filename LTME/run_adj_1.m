function error = run_adj

%
% Aug. 2007, Y.H. Kim
% For the summer school, solution of parameter estimation.
% based on secant method
%

error_ratio = 1.2;

fom.p = 10.0*error_ratio;
fom.r = 32.0;
fom.b = 2.66666667;

alpha = 2*10^-6;

dt = 0.01;
nstop = 200;

load obs.mat;

fom.x(1) = 1*error_ratio;
fom.y(1) = 3*error_ratio;
fom.z(1) = 5*error_ratio;
adj.x(1) = 0;
adj.y(1) = 0;
adj.z(1) = 0;
adj.p(1) = 0;
adj.r(1) = 0;
adj.b(1) = 0;

for iter=1:200
        
    fom.x(1) = fom.x(1) - alpha*adj.x(1);
    fom.y(1) = fom.y(1) - alpha*adj.y(1);
    fom.z(1) = fom.z(1) - alpha*adj.z(1);
    fom.p = fom.p - alpha*adj.p(1);
    fom.r = fom.r - alpha*adj.r(1);
    fom.b = fom.b - alpha*adj.b(1);

    for i=1:nstop
      [xp,yp,zp] = florenz(fom.x(i),fom.y(i),fom.z(i),fom.p,fom.r,fom.b);
      fom.x(i+1) = fom.x(i)+dt*xp;
      fom.y(i+1) = fom.y(i)+dt*yp;
  	  fom.z(i+1) = fom.z(i)+dt*zp;
    end
    
    adj.x(nstop+1) = 0;
    adj.y(nstop+1) = 0;
    adj.z(nstop+1) = 0;
    adj.p(nstop+1) = 0;
    adj.r(nstop+1) = 0;
    adj.b(nstop+1) = 0;

    [error,obs_index] = size(obs.pos);

    for i=nstop:-1:1
    
        if((i+1) == obs.pos(obs_index))
            innovx(obs_index) = -obs.x(obs_index) + fom.x(i+1);
            innovy(obs_index) = -obs.y(obs_index) + fom.y(i+1);
            innovz(obs_index) = -obs.z(obs_index) + fom.z(i+1);
            adj.x(i+1) = adj.x(i+1) + 10*(fom.x(i+1) - obs.x(obs_index));
            adj.y(i+1) = adj.y(i+1) + 10*(fom.y(i+1) - obs.y(obs_index));
            adj.z(i+1) = adj.z(i+1) + 10*(fom.z(i+1) - obs.z(obs_index));
            obs_index = max(obs_index - 1,1);
        end
    
      [xa,ya,za,pe,re,be] = alorenz(adj.x(i+1),adj.y(i+1),adj.z(i+1),...
          fom.x(i),fom.y(i),fom.z(i),fom.p,fom.r,fom.b);
       adj.x(i) = adj.x(i+1)+dt*xa;
       adj.y(i) = adj.y(i+1)+dt*ya;
       adj.z(i) = adj.z(i+1)+dt*za;
       adj.p(i) = adj.p(i+1)+dt*pe;
       adj.r(i) = adj.r(i+1)+dt*re;
       adj.b(i) = adj.b(i+1)+dt*be;
       
    end

    costf(iter) = innovx*innovx'+innovy*innovy'+innovz*innovz';
    fom_ana(iter) = fom;
end
clf;
plot(costf);
save fom_ana.mat fom_ana;
    