function f = dtlz2a(x)

alpha = 1;

%% g function
sum = 0;
for i = 3:8
sum = sum + (x(:,i)-0.5)^2;
end
g = sum;

%% cost functions
f(:,1) = (1+g)*cos(((x(:,1)^alpha)*pi)/2)*cos(((x(:,2)^alpha)*pi)/2);
f(:,2) = (1+g)*cos(((x(:,1)^alpha)*pi)/2)*sin(((x(:,2)^alpha)*pi)/2);
f(:,3) = (1+g)*sin(((x(:,1)^alpha)*pi)/2);
end
