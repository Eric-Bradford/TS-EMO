function f = vlmop2(x)

dim = 2;

transl = 1./sqrt(dim);
part1 = (x(:,1) - transl).^2 + (x(:,2) - transl).^2;
part2 = (x(:,1) + transl).^2 + (x(:,2) + transl).^2;

f(:,1) = 1 - exp( -part1 );
f(:,2) = 1 - exp( -part2 );
