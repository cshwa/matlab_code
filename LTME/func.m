function F = func(x)

% function F = func(x)
% Initial state estimation problem
%
% Run the Lorenz Model and calculate the constfunction
% Programmed by Y.H. Kim for the summer school on Aug. 2007 
%
% x : input state vector
% F : output cost function
%

p = 10.0;
r = 32.0;
b = 2.66666667;

dt = 0.01;
nstop = 200;

fom.x(1) = x(1);
fom.y(1) = x(2);
fom.z(1) = x(3);

for i=1:nstop
    [xp,yp,zp] = florenz(fom.x(i),fom.y(i),fom.z(i),p,r,b);
    fom.x(i+1) = fom.x(i)+dt*xp;
    fom.y(i+1) = fom.y(i)+dt*yp;
    fom.z(i+1) = fom.z(i)+dt*zp;
end

load obs.mat;
innovx = fom.x(obs.pos) - obs.x;
innovy = fom.y(obs.pos) - obs.y;
innovz = fom.z(obs.pos) - obs.z;

F = innovx*innovx'+innovy*innovy'+innovz*innovz';
return;