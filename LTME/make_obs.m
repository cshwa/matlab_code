function error = make_obs

% function error = make_obs
% make observation
% Programmed by Y.H. Kim for the summer school on Aug. 2007 
%

p = 10.0;
r = 32.0;
b = 2.66666667;

x = 1;
y = 3;
z = 5;

dt = 0.01;
nstop = 200;

fom.x(1) = x;
fom.y(1) = y;
fom.z(1) = z;

obs_index = 0;

rx = 0*random('normal',0,0.1,1,100); % observation error
ry = 0*random('normal',0,0.1,1,100); % observation error
rz = 0*random('normal',0,0.1,1,100); % observation error

for i=1:nstop
    [xp,yp,zp] = florenz(fom.x(i),fom.y(i),fom.z(i),p,r,b);
    fom.x(i+1) = fom.x(i)+dt*xp;
    fom.y(i+1) = fom.y(i)+dt*yp;
    fom.z(i+1) = fom.z(i)+dt*zp;
    if(mod(i,5) == 0)
        obs_index = obs_index + 1;
        obs.pos(obs_index) = i+1;
        obs.x(obs_index) = fom.x(i+1)+rx(obs_index);
        obs.y(obs_index) = fom.y(i+1)+ry(obs_index);
        obs.z(obs_index) = fom.z(i+1)+rz(obs_index);
    end
end

save obs.mat obs;

    