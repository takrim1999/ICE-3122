x = [-2*pi:0.1:2*pi];
y1=sin(x);
y2=cos(x);
y3=tan(x);
y4=sec(x);
y5=csc(x);
y6=cot(x);
axis([-2*pi 2*pi -2*pi 2*pi]);
hold on
plot(x,y1);
plot(x,y2);
plot(x,y3);
plot(x,y4);
plot(x,y5);
plot(x,y6);
legend('sin','cos','tan','sec','cosec','cot');
hold off
