function error = run_adj

p = 10.0;
r = 32.0;
b = 2.66666667;

error_ratio = 1.3;
alpha = 2*10^-5;

dt = 0.01;
nstop = 200;

load obs.mat;

fom.x(1) = 1*error_ratio;
fom.y(1) = 3*error_ratio;
fom.z(1) = 5*error_ratio;
adj.x(1) = 0;
adj.y(1) = 0;
adj.z(1) = 0;

for iter=1:100
        
    fom.x(1) = fom.x(1) - alpha*adj.x(1);
    fom.y(1) = fom.y(1) - alpha*adj.y(1);
    fom.z(1) = fom.z(1) - alpha*adj.z(1);

    for i=1:nstop
      [xp,yp,zp] = florenz(fom.x(i),fom.y(i),fom.z(i),p,r,b);
      fom.x(i+1) = fom.x(i)+dt*xp;
      fom.y(i+1) = fom.y(i)+dt*yp;
  	  fom.z(i+1) = fom.z(i)+dt*zp;
    end
    
    adj.x(nstop+1) = 0;
    adj.y(nstop+1) = 0;
    adj.z(nstop+1) = 0;

    [error,obs_index] = size(obs.pos);

    for i=nstop:-1:1
    
        if((i+1) == obs.pos(obs_index))
            innov(obs_index) = -obs.x(obs_index) + fom.x(i+1);
            adj.x(i+1) = adj.x(i+1) + 10*(fom.x(i+1) - obs.x(obs_index));
            adj.y(i+1) = adj.y(i+1) + 10*(fom.y(i+1) - obs.y(obs_index));
            adj.z(i+1) = adj.z(i+1) + 10*(fom.z(i+1) - obs.z(obs_index));
            obs_index = max(obs_index - 1,1);
        end
    
      [xa,ya,za] = alorenz(adj.x(i+1),adj.y(i+1),adj.z(i+1),fom.x(i),fom.y(i),fom.z(i),p,r,b);
       adj.x(i) = adj.x(i+1)+dt*xa;
       adj.y(i) = adj.y(i+1)+dt*ya;
       adj.z(i) = adj.z(i+1)+dt*za;
    end

    costf(iter) = innov*innov';
    fom_ana(iter) = fom;
    disp([sprintf('%10.4e %10.4e %10.4e',fom.x(1),fom.y(1),fom.z(1))])
end
clf;
plot(costf);
set(gcf,'position',[0 0 500 400]);
xlabel('assimilation step','fontsize',15);
ylabel('cost function','fontsize',15);
save fom_ana.mat fom_ana;
    
    