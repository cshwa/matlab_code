function F_prime = dfunc(x)

% function F_prime = dfunc(x)
% Parameter estimation problem
%
% Run the adjoint Lorenz Model and calculate the gradient of const function
% Programmed by Y.H. Kim for the summer school on Aug. 2007 
%
% x : input state vector
% F_prime : output gradient of cost function
%

dt = 0.01;
nstop = 200;

load obs.mat;

fom.x(1) = x(1);
fom.y(1) = x(2);
fom.z(1) = x(3);
p = x(4);
r = x(5);
b = x(6);

for i=1:nstop
  [xp,yp,zp] = florenz(fom.x(i),fom.y(i),fom.z(i),p,r,b);
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
        adj.x(i+1) = adj.x(i+1) + (fom.x(i+1) - obs.x(obs_index));
        adj.y(i+1) = adj.y(i+1) + (fom.y(i+1) - obs.y(obs_index));
        adj.z(i+1) = adj.z(i+1) + (fom.z(i+1) - obs.z(obs_index));
        obs_index = max(obs_index - 1,1);
    end
    
   [xa,ya,za,pe,re,be] = alorenz(adj.x(i+1),adj.y(i+1),adj.z(i+1),fom.x(i),...
       fom.y(i),fom.z(i),p,r,b);
   adj.x(i) = adj.x(i+1)+dt*xa;
   adj.y(i) = adj.y(i+1)+dt*ya;
   adj.z(i) = adj.z(i+1)+dt*za;
   adj.p(i) = adj.p(i+1)+dt*pe;
   adj.r(i) = adj.r(i+1)+dt*re;
   adj.b(i) = adj.b(i+1)+dt*be;
end

F_prime = [adj.x(1) adj.y(1) adj.z(1) adj.p(1) adj.r(1) adj.b(1)];

save dfunc.mat x F_prime;

return;
    