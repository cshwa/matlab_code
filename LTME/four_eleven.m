clear all
i=2
x=pi/4
true=cos(pi/4)
epsilons=0.5

app(1)=1
tre(1)=(true-app(1))/true *100
are=1
aree(1)=1
while epsilons<are;
app((i/2)+1)= (app(i/2) + ((-1)^(i/2))*(x^(i))/factorial(2^(i/2)))

are=abs((app((i/2)+1)-app(i/2))/app((i/2)+1)) *100
aree((i/2)+1)=are
tre((i/2)+1)=(true-app((i/2)+1))/true * 100

i=i+2
end