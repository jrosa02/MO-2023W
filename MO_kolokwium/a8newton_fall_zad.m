clear all
close all
syms x1 x2
f = x1-(4*(x1+x2).*exp(-x1.^2 - x2.^2)+0.1)^2;
f_num = matlabFunction(f);

figure;
[X, Y] = meshgrid(-3:0.1:3, -3:0.1:3);
Z = f_num(X, Y);
contour(X, Y, Z, 30);

P01=[1, 0.5];
eps = 1e-3;
Nmax = 1000;

[x, i, x_vec] = max_fall(f, P01, eps, Nmax)
plot_result(f_num, x, x_vec)
disp([x, i, x_vec])

[x, i, x_vec] = newton_fall(f, P01, eps, Nmax)
plot_result(f_num, x, x_vec)
disp([x, i, x_vec])

%% funckje właściwe

function [x, i, x_vec] = max_fall(f, x0, eps, Nmax)
i = 1;
x_vec = [x0];
x = x0;
x_prev = x0;
syms x1 x2
grad = gradient(f, [x1, x2]).';
f_num = matlabFunction(f);

while i <= Nmax
    x_prev = x;
    d = -double(subs(grad, [x1, x2], [x(1), x(2)]));
    f_l = @(h) f_num(x(1) + h*d(1), x(2) + h*d(2));
    h = fminsearch(f_l, 0);
    x = x_prev + h*d;
    x_vec(i, :) = x;
    i = i + 1;
    if rms(x_prev - x) <= eps
        return;
    end %if
end %while
end %function

function [x, i, x_vec] = newton_fall(f, x0, eps, Nmax)
i = 1;
x_vec = [x0];
x = x0;
x_prev = x0;
syms x1 x2
grad = gradient(f, [x1, x2]).';
hesjan = hessian(f, [x1, x2]);
inv_hesjan = inv(hesjan);
f_num = matlabFunction(f);

while i <= Nmax
    x_prev = x;
    f_n = -inv_hesjan*grad.';
    d = -double(subs(f_n, [x1, x2], [x(1), x(2)]));
    f_l = @(h) f_num(x(1) + h*d(1), x(2) + h*d(2));
    h = fminsearch(f_l, 0);
    x = x_prev + h*d.';
    x_vec(i, :) = x;
    i = i + 1;
    if rms(x_prev - x) <= eps
        return;
    end %if
end %while
end %function

%% wizualizacja

function [] = plot_result(f, x, x_vec)
figure;
[X, Y] = meshgrid(-3:0.1:3, -3:0.1:3);
Z = f(X, Y);
contour(X, Y, Z, 100);
title('Contour Plot of the Objective Function');
xlabel('x1');
ylabel('x2');
hold on;

scatter(x_vec(:, 1), x_vec(:, 2), 10, 'r', 'filled');
hold on; % to keep the scatter plot
line(x_vec(:, 1), x_vec(:, 2));
hold off;
text(x(1), x(2), '  Koniec', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

legend('Funkcja Celu', 'Punkt Końcowy');
grid on;
hold off;

end
