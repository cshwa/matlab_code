t = [1:1000];
a = sin(t/20)*rand(1,1000);
a = sin(t/20)+rand(1,1000);
b = cos(t/20)*rand(1,1000);
b = cos(t/20)+rand(1,1000);
wtc(a,b)
help wtc
wtc(a,b,'MonteCarloCount',1)
edit wtc
wtc(a,b,'MonteCarloCount',1)
b = cos(t/20+pi/2)+rand(1,1000);
plot(a)
hold on
plot(b)
wtc(a,b,'MonteCarloCount',1)
b = cos(t/20-pi/2)+rand(1,1000);
wtc(a,b,'MonteCarloCount',1)
b = cos(t/20)+rand(1,1000);
wtc(a,b,'MonteCarloCount',1)
a = sin(t/20-pi/2)+rand(1,1000);
wtc(a,b,'MonteCarloCount',1)
a = sin(t/20+pi/2)+rand(1,1000);
wtc(a,b,'MonteCarloCount',1)
b = cos(t/20+pi/2)+rand(1,1000);
wtc(a,b,'MonteCarloCount',1)