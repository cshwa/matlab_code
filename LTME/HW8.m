clear all
c=1;
dt=0.08;
dx=1;
u(2:10,1)=0;
u(26:35,1)=0;
for x=10:18
    u(x,1)=x/8-10/8;
end
for x=18:26
    u(x,1)=-x/8+26/8;
end
u0=u;
u2=u;
for x=1:33
    f(x)=1/x;
end
for t=1:32
    for x=1:33
        if (x-t>0)
            u(x,t+1)=u(x-t,1);
        else
            u(x,t+1)=u(x-t+33,1);
        end
    end
end

for t=1:33
    plot(u(:,t))
    ylim([0 1]);
    xlim([0 35]);
    pause(0.1)
    
end

% 2nd order centered space differecting
t=1
for x=2:34
%         u2(x,t+1)=dt/2/dx*(u2(x+1,t) -2*u2(x,t) + u2(x-1,t)) +u2(x,t);
%         if (x==2)
%             u2(x,t+1)=dt/dx*(u2(3,t) -u2(34,t)) +u2(2,t-1);
%         elseif(x==34)
%             u2(x,t+1)=dt/dx*(u2(2,t) -u2(33,t)) +u2(34,t-1);
%         else
%             u2(x,t+1)=dt/dx*(u2(x+1,t) -u2(x-1,t)) +u2(x,t-1);
%         end     
        if (x==2)
            u2(x,t+1)=-dt/dx*(u2(x+1,t) -u2(2,t)) +u2(2,t);
        elseif(x==34)
            u2(x,t+1)=-dt/dx*(u2(2,t) -u2(34,t)) +u2(34,t);
        else
            u2(x,t+1)=-dt/dx*(u2(x+1,t) -u2(x,t)) +u2(x,t);
        end
    end
for t=2:400
    for x=2:34
%         u2(x,t+1)=dt/2/dx*(u2(x+1,t) -2*u2(x,t) + u2(x-1,t)) +u2(x,t);
        if (x==2)
            u2(x,t+1)=-dt/dx*(u2(3,t) -u2(34,t)) +u2(2,t-1);
        elseif(x==34)
            u2(x,t+1)=-dt/dx*(u2(2,t) -u2(33,t)) +u2(34,t-1);
        else
            u2(x,t+1)=-dt/dx*(u2(x+1,t) -u2(x-1,t)) +u2(x,t-1);
        end     
%         if (x==2)
%             u2(x,t+1)=dt/dx*(u2(x+1,t) -u2(2,t)) +u2(2,t);
%         elseif(x==34)
%             u2(x,t+1)=dt/dx*(u2(2,t) -u2(34,t)) +u2(34,t);
%         else
%             u2(x,t+1)=dt/dx*(u2(x+1,t) -u2(x,t)) +u2(x,t);
%         end
    end
end

% for t=1:40
%     plot(u2(:,t))
%     ylim([0 1]);
%     xlim([1 34]);
%     pause(0.1)
% end
hold on
plot(u2(:,400))
ylim([0 1]);
xlim([1 34]);



 %FE Method
  %A matrix
  A = zeros(33);
  for i = 1:33
      if i == 1
          A(i,i) = 4;
          A(i,2) =1;
          A(i,33) = 1; %%%%%%%%%%%%
      elseif i ==33
          A(i,i) = 4;
          A(i,1) = 1; %%%%%%%%%%%%
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
            B(i,33) =-1;   %%%%%%%
        elseif i==33
            B(i,1) = 1;   %%%%%%%%%%
            B(i,32) = -1;
        else
            B(i,i+1) = 1;
            B(i,i-1) = -1;
        end
   end
   AB = inv_A * B;
  %cal 
u_fem(:,1)=u0(2:34,:);
u_fem(:,2)=u2(2:34,2);
for i = 2 : 400
       u_fem(:,i+1) = u_fem(:,i-1) - 6*dt * AB * u_fem(:,i); 
end
plot(u_fem(:,400))
ylim([0 1]);
xlim([1 34]);

%3) spectral method

 %make initial condition
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

 u_spe = zeros(33,400) ;
 u_spe(:,1)=u0(2:34,:);
u_spe(:,2)=u2(2:34,2);

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
 u_spe(:,1)=u0(2:34,:);
u_spe(:,2)=u2(2:34,2);

% figure;
plot(u_spe(:,400)+0.25)
ylim([0 1]);
xlim([1 34]);

u22 = sum(u(:,1).*u(:,1))/2;
u222 = sum(u2(:,1).*u2(:,1))/2;
u_fem2 = sum(u_fem(:,1).*u_fem(:,1))/2;
u_spe2 = sum(u_spe(:,1).*u_spe(:,1))/2;