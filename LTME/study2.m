clear all

x=1:1000
t(1:1000)=1000
A0=1
lambda1=100;
T1=60;

lambda2=100;
T2=50;
for i=1:1000
    Scenario1(i)=A0 * ( sin( 2*pi *(x(i)/lambda1 - t(i)/T1) ) + sin( 2*pi *(x(i)/lambda2 - t(i)/T2) ) );
end
subplot(3,2,1)
plot(x,Scenario1)
title('Scenario1')

lambda2=90;
T2=60;
for i=1:1000
    Scenario2(i)=A0 * ( sin( 2*pi *(x(i)/lambda1 - t(i)/T1) ) + sin( 2*pi *(x(i)/lambda2 - t(i)/T2) ) );
end
subplot(3,2,2)
plot(x,Scenario2)
title('Scenario2')

lambda2=90;
T2=50;
for i=1:1000
    Scenario3(i)=A0 * ( sin( 2*pi *(x(i)/lambda1 - t(i)/T1) ) + sin( 2*pi *(x(i)/lambda2 - t(i)/T2) ) );
end
subplot(3,2,3)
plot(x,Scenario3)
title('Scenario3')

lambda2=100;
T2=-60;
for i=1:1000
    Scenario4(i)=A0 * ( sin( 2*pi *(x(i)/lambda1 - t(i)/T1) ) + sin( 2*pi *(x(i)/lambda2 - t(i)/T2) ) );
end
subplot(3,2,4)
plot(x,Scenario4)
title('Scenario4')

lambda2=50;
T2=-30;
for i=1:1000
    Scenario5(i)=A0 * ( sin( 2*pi *(x(i)/lambda1 - t(i)/T1) ) + sin( 2*pi *(x(i)/lambda2 - t(i)/T2) ) );
end
subplot(3,2,5)
plot(x,Scenario5)
title('Scenario5')

lambda2=95;
T2=-30;
for i=1:1000
    Scenario6(i)=A0 * ( sin( 2*pi *(x(i)/lambda1 - t(i)/T1) ) + sin( 2*pi *(x(i)/lambda2 - t(i)/T2) ) );
end
subplot(3,2,6)
plot(x,Scenario6)
title('Scenario6')