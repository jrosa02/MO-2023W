%% przykładowe gówno
f = @(x) 0.01*sin(10*x(1))+(0.07*x(1).^2-0.4).^2;
ft = @(x) f(x + 32);

px = [-10 -5 0 5 10];
py = [0.0001 -0.1 0 0.1 -0.0001];
stopien=5
W = polyfit(px, py, stopien)
fx=@(x)(0.05*sin(1*x)+W(1)*x.^stopien+W(2)*x.^(stopien-1)+W(3)*x.^(stopien-2)+W(4)*x.^(stopien-3)+W(5)*x.^(stopien-4)+W(6)*x.^(stopien-5))*(0.001*sin(x/10))+0.000001*x

x0 = -0.5;
zad2(f, x0, x0+0.1, 1.5, 0.5, 100)

%Warunek stopu jest konieczny do dziania algorytmów numerycznych
%Zatrzymuje działanie algorytmu gdy zaczyna oddalać się od rozwiązania
%Algorytm kontrakcji jest dobrym uzupełnieniem algorytmu ekspansji

f = fx;
x0 = -0.5;
zad2(f, x0, x0+0.1, 1.5, 0.5, 100)

%% funkcje właściwe
function [a, b, i] = contraction(f, beta, x0, x1, Nmax)
if f(x0) > f(x1)
    xtemp = x0;
    x0 = x1;
    x1 = xtemp;
end

i = 1;
xj = x1;

while f(x0) < f(xj)
    if i > Nmax
        a = -1.*inf;
        b = inf;
        return;
    end
    xj = x0 + beta.^(i)*(x1 - x0);
    i = i +1;
end
a = x0;
b = xj;
return
end


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
%% funckje wyświetlające

function [] = zad2(f, x0, x1, alpha, beta, Nmax)
disp("----------------START-------------------")
f
x0
x1
alpha
beta
[a_exp, b_exp, i_exp] = expansion(f, x0, x1, alpha, Nmax)
[a_contr, b_contr, i_contr] = contraction(f, beta, a_exp, b_exp, Nmax)
xmin = min([x0, x1, a_exp, b_exp, a_contr, b_contr]);
xmax = max([x0, x1, a_exp, b_exp, a_contr, b_contr]);
xlimits = [xmin, xmax];
figure;
hold on;
fplot(f, xlimits);
plot([x0, x1], [f(x0), f(x1)], '|c',...
    [a_exp, b_exp], [f(a_exp), f(b_exp)], '*r',...
    [a_contr, b_contr], [f(a_contr), f(b_contr)], '+b');
legend(["Func", "Start", "Exp", "Contr"]);
hold off;
disp("----------------END---------------------")
end


function [] = zad2contr(f, x0, x1, beta, Nmax)
disp("----------------START-------------------")
f
x0
x1
beta
[a_contr, b_contr, i_contr] = contraction(f, beta, x0, x1, Nmax)
xmin = min([x0, x1, a_contr, b_contr]);
xmax = max([x0, x1, a_contr, b_contr]);
xlimits = [xmin, xmax];
figure;
hold on;
fplot(f, xlimits);
plot([x0, x1], [f(x0), f(x1)], '|c',...
    [a_contr, b_contr], [f(a_contr), f(b_contr)], '+r');
%legend(["Func", "Start", "Contr"]);
hold off;
disp("----------------END---------------------")
end