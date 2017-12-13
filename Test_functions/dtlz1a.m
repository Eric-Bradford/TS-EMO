function f = dtlz1a(x)

%% g function
sum = 0;
for i = 2:6
sum = sum + (x(:,i)-0.5)^2 - cos(2*pi*(x(:,i)-0.5));
end
g = 100*(5+sum);

%% cost functions
f(:,1) = 0.5.*x(:,1)*(1+g);
f(:,2) = 0.5*(1-x(:,1))*(1+g);
end
