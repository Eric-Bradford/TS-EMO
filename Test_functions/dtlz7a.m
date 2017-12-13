function f = dtlz7a(x)

%% cost functions
f(:,1) = x(:,1);
f(:,2) = x(:,2);

%% g function
sum1 = 0;
sum2 = 0;
for i = 3:8
sum1 = sum1 + x(:,i);
end
g = 1+(9/6)*sum1;
for i = 1:2
sum2 = sum2 + (f(:,i)/(1+g))*(1+sin(3*pi*f(:,i)));
end
h = 3 - sum2;

%% cost functions
f(:,3) = (1+g)*h;

