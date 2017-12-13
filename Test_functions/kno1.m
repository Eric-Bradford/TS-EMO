function f = kno1(x)

r = 9 - (  3*sin( 5/2*(x(1)+ x(2))^2 ) + 3*sin( 4*(x(1) + x(2)) ) + 5*sin( 2*(x(1) + x(2)) + 2 ) );

phi = pi/12 * (x(1) - x(2) + 3) ;

f(:,1) = 20 - r * cos(phi) ;
f(:,2) = 20 - r * sin(phi);
