% arrow.m È°¿ë
% m-file: colorarrow_exam.m
x = 0:10;
y = sin(x);
plot(x,y,'o-')
axis([0 10 -4 2])
for n = 1:length(x)
 arrow([x(n) -4],[x(n) y(n)],10+4*n,...
 10+5*n,5+3*n,'EdgeColor','r','FaceColor','r')
end
title('color arrow example')