function[a, b]= fit_power_law(x, y)
if(x(1)<0.5)
    x = x + 0.5;
end;
y(y<1e-5) = 1e-5;

n = length(x);
[~, ii] = min((x/max(x)).^2+(y/max(y)).^2);

x = log(x(1:ii));
y = log(y(1:ii));
c = polyfit(x, y, 1);
a = exp(c(1));
b = c(2);