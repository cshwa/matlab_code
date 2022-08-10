clc;clear all;close all;
%HW #8
clc;clear all;close all;
dt = 0.1 ;
dx = 1.0 ;
%1) FD Method
 %set values
u_fdm = zeros(33,400) ; %u = u(x,t)
u_fem = zeros(33,400) ;
u_spe = zeros(33,400) ;

 %set initial value
 u_ini = zeros(33,2);
 for i = 1 : 8
     u_ini(i+9,1) = 1/8*i ;
     u_ini(i+17,1) = 1-1/8*i ;
 end
 
 %FDM of forward differencing for initial value
 for i = 1 :33
     if i ==1
         u_ini(i,2) = u_ini(i,1) - dt/dx/2*( u_ini(i,1)-u_ini(32,1)) ;
     elseif i==33
         u_ini(i,2) = u_ini(i,1) - dt/dx/2*( u_ini(i,1) - u_ini(32,1)) ;
     else
        u_ini(i,2) = u_ini(i,1) - dt/dx/2*(u_ini(i,1)-u_ini(i-1,1));
     end
end
 
 %apply initial value
 u_fdm(:,1:2) = u_ini;
 u_fem(:,1:2) = u_ini;
 u_spe(:,1:2) = u_ini;
 
 for n = 2:399
     for m = 1:33
        if m ==1
            u_fdm(m,n+1) = u_fdm(m,n-1) - dt/dx * (u_fdm(m+1,n) -u_fdm(32,n)); 
        elseif m==33
            u_fdm(m,n+1) = u_fdm(m,n-1) - dt/dx * (u_fdm(2,n) -u_fdm(m-1,n));
        else
            u_fdm(m,n+1) = u_fdm(m,n-1) - dt/dx * (u_fdm(m+1,n) -u_fdm(m-1,n));
        end
     end
 end
 
 %FE Method
  %A matrix
  A = zeros(33);
  for i = 1:33
      if i == 1
          A(i,i) = 4;
          A(i,2) =1;
          A(i,32) = 1;
      elseif i ==33
          A(i,i) = 4;
          A(i,2) = 1;
          A(i,32) = 1;
      else
          A(i,i) = 4;
          A(i,i-1) = 1;
          A(i,i+1) = 1;
      end
  end
  inv_A = inv(A);
   %B matrix
   B = zeros(33);
   for i = 1:33
        if i ==1
            B(i,i+1) = 1;
            B(i,32) =-1;
        elseif i==33
            B(i,2) = 1;
            B(i,32) = -1;
        else
            B(i,i+1) = 1;
            B(i,i-1) = -1;
        end
   end
   AB = inv_A * B;

   %calculate
for i = 2 : 399
       u_fem(:,i+1) = u_fem(:,i-1) - 6*dt * AB * u_fem(:,i); 
end
   
%3) spectral method
C = zeros(33,33);D = zeros(33,400);
for i = 1 : 33
    b = [];c=[];
    for j = 1 : 16
        a = [1];
        b = [b cos(2*j*pi/32* (i-1))];
        c = [c sin(2*j*pi/32 * (i-1))];
    end
    C(i,:) = [a b c] ;
end

 %make initial condition
for i = 1:2 
    D(:,i) = inv(C) * u_spe(:,i);
end

for i = 2: 399
    D(1,i) =D(1); 
    for j = 1 : 16
        x = j+1 ;
        D(x,i+1) = (2-(dt)^2 *(2*j*pi/32)^2)*D(x,i) -D(x,i-1);
        x = j+17;
        D(x,i+1) = (2-(dt)^2 *(2*j*pi/32)^2)*D(x,i) -D(x,i-1);
    end
end

u_spe = C * D ;

%plot 
x = [0:32];
plot(x,u_ini(:,1),'k','Linewidth',2);hold on;
plot(x,u_fdm(:,320),'-.r','Linewidth',1.5);
plot(x,u_fem(:,320),'-.g','Linewidth',1.5);
plot(x,u_spe(:,320),'-.b','Linewidth',2);
legend('Analysis','Finite difference','Finite element','Spectral', 'location', 'northwest');
set(gca,'xtick', [0:4:32]); xlim([0 32])
saveas(gcf,'hw8.png','png');

%u^2 
u_ini2 = sum(u_ini(:,1).*u_ini(:,1))/2;
u_fdm2 = sum(u_fdm(:,1).*u_fdm(:,1))/2;
u_fem2 = sum(u_fem(:,1).*u_fem(:,1))/2;
u_spe2 = sum(u_spe(:,1).*u_spe(:,1))/2;