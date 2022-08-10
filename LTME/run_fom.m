function error = run_fom(x)

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

plot(fom.y);
axis([0 200 -20 40]);
xlabel('time step','fontsize',15);
ylabel('x','fontsize',15);
set(gcf,'position',[0 0 800 400]);
save fom.mat fom;
