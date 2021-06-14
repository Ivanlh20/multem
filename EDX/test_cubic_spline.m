close all
clear
clf reset

ff = @(x) sin(x);
x = [-2 0 2 4 6 7 8 9 10 11 12];
for i = 1: length(x)
    y(i) = ff(x(i));
end

n = length(x);
h = zeros(n - 1, 1);
for i = 1: n - 1
    h(i) = x(i + 1) - x(i);
end

a = zeros(n, 1);
b = zeros(n, 1);
d = zeros(n, 1);
f = zeros(n, 1);

for j = 2: n - 1
    a(j) = 2*(h(j) + h(j - 1));
    b(j) = h(j - 1);
    d(j) = h(j);
    f(j) = 6*(y(j + 1) - y(j)/h(j)) - (6 *(y(j) - y(j - 1)) / h(j - 1));
end

% Starting conditions
a(1) = 3;
f(1) = 0;

a(n) = 0;
f(n) = 0;

% Coefficents
c = tridiag(a, b, d, f);

t = linspace (x(1), x(n), 301);
for k = 1: length(t)
    tk = t(k);
    y_s(k) = spline_aux(x, y, c, tk);
    z(k) = ff(tk);
end

y_m = spline(x, y, t);
figure(1); clf;
plot(t, z, 'b-', 'linewidth', 2)
hold on;
plot(t, y_s, 'r.', t, y_m, 'k.', x, y, 'go');

function [x] = tridiag(a, b, d, f)
    n = length(f);
    alfa = zeros(n, 1);
    beta = zeros(n, 1);

    alfa(1) = a(1);
    for i = 2: n
        beta(i) = b(i) / alfa(i - 1);
        alfa(i) = a(i) - (beta(i)*d(i - 1));
    end

    y(1) = f(1);
    for i = 2:n
        y(i) = f(i) - beta(i)*y(i - 1);
    end

    x(n) = y(n) / alfa(n);
    for i = (n-1):-1:1
        x(i) = (y(i) - (d(i)*x(i + 1))) / alfa(i);
    end
end

function [s] = spline_aux(x, w, c, tk)
    n = length(x);

    h = zeros(n - 1, 1);
    for i = 1: n - 1
        h(i) = x(i+1) - x(i);
    end

     for i = 1: n - 1
         if (x(i) <= tk && tk <= x(i+1))
            break
         end
     end

    s1 = c(i)*((x(i+1) - tk)^3)/(6*h(i));
    s2 = c(i+1)*((tk - x(i))^3)/(6*h(i));
    s3 = (w(i)/h(i) - (c(i)*h(i)/6))*(x(i+1) - tk);
    s4 = (w(i+1)/h(i) - (c(i+1)*h(i)/6))*(tk - x(i));
    s = s1 + s2 + s3 + s4;
end