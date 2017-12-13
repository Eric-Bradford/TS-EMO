function f = oka2( x )

f(:,1) = x(1);
f(:,2) = 1 - 1/(4*pi^2) * (x(1)+pi)^2 + (abs( x(2) - 5*cos(x(1)) ))^(1/3) + (abs( x(3) - 5*sin(x(1)) ))^(1/3) ;