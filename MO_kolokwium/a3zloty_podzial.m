%Funkcja 1
f= @(x) 0.01*sin(10*x(1))+(0.07*x(1).^2-0.4).^2

x0 = -0.5;
zad3(f, x0, x0+0.1, 2, epsilon, 1000)

%% wizuaizacja
function [] = zad3(f, x0, x1, alpha, eps, Nmax)
disp("----------------START-------------------")
f
x0
x1
eps
a = x0; b = x1;
[x, i] = goldendiv(f, a, b, eps, Nmax)
xmin = min([x0, x1, a, b]);
xmax = max([x0, x1, a, b]);
xlimits = [xmin, xmax];
figure;
hold on;
fplot(f, xlimits);
plot([x0, x1], [f(x0), f(x1)], '|c',...
       [a, b], [f(a), f(b)], '*r',...
    [x], [f(x)], '+k');
legend(["Func", "Start", "Exp", "GoldDiv"]);
hold off;
disp("----------------END---------------------")
end

%% funkcje właściwe
function [x, i] = goldendiv(f, a, b, eps, Nmax)
if a > b
    temp = a;
    a = b;
    b = temp;
end

i = 0;

if a == b
    return;
end

i = 1;

alpha = 0.6180339;
a_v = [a];
b_v = [b];
c_v = [b - alpha*(b - a)];
d_v = [a + alpha*(b - a)];

while ( abs(b_v(i) - a_v(i)) ) >= eps
    if f(c_v(i)) < f(d_v(i))
        a_v(i+1) = a_v(i);
        b_v(i+1) = d_v(i);
        c_v(i+1) = b_v(i+1) - alpha*(b_v(i+1) - a_v(i+1));
        d_v(i+1) = c_v(i);
    else
        a_v(i+1) = c_v(i);
        b_v(i+1) = b_v(i);
        d_v(i+1) = a_v(i+1) + alpha*(b_v(i+1) - a_v(i+1));
        c_v(i+1) = d_v(i);
    end
    i = i + 1;
    
    if i > Nmax
        x = inf;
        return;
    end
end
    i = i -1;
    x = ( a_v(i) + b_v(i) ) /2;
end


%% expansja
function [a, b, i] = expansion(f, xzero, x1, alpha, Nmax)
x0 = xzero;
i = 0;
if f(x0) == f(x1)
    a = x0;
    b = x1;
    return
end

if f(x1) > f(x0)
    x1 = -1*(x1 - xzero) + xzero;
    if f(x1) == f(-1*(x1 - xzero) + xzero)
        a = -1*(x1 - xzero) + xzero;
        b = x1;
        return
    end
end

x2 = x1;
x1 = x0;
while f(x1) > f(x2)
    if i > Nmax
        a = -inf;
        b = inf;
        return 
    end
    i = i + 1;
    x0 = x1;
    x1 = x2;
    x2 = alpha * (x2 - xzero) + xzero;
end

if x0 < x2
    a = x0;
    b = x2;
    return
end
a = x2;
b = x0;
return;
end