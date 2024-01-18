f= @(x1, x2) 0.3*x1+0.1*x2+(-3.5+0.5*x1.^2+0.5*x2.^2).^2+100*x1 .* exp(-x1.^2 - x2.^2)
P0=[1.5, 0];
P1=[-2,-2];
e_n1 = [1 -1; 1 1];
[x, it] = Powell(f, P0, e_n, n, Nmax, epsilon)
plot_result(f, x)

%% Problemowe 

global G1;
global t;
s = tf([1 0], [0 1]);

G1=tf([0 1],[1 2 2.25 1.25]);
[G1p_Licz,G1p_Mian]=pade(0.5,3);
G1p=tf(G1p_Licz,G1p_Mian);
G1=G1*G1p
step(G1)

[y, t] = step(G1, t);

t = linspace(0, 25, 5000);

param0 = [6, 0];

out_par = fminsearch(@f_obj, param0)

G_approx = tf([0, out_par(2)], [out_par(1), 1]);


plot(t, step(G1, t), '-r', t, step(G_approx, t), '--b')

%% wizualizacja

function [] = plot_result(f, x)
figure;
[X, Y] = meshgrid(-3:0.1:3, -3:0.1:3);
Z = f(X, Y);
contour(X, Y, Z, 50);
title('Contour Plot of the Objective Function');
xlabel('x1');
ylabel('x2');
hold on;

% Oznacz punkt końcowy na wykresie konturów
scatter(x(end, 1), x(end, 2), 100, 'r', 'filled');
text(x(end, 1), x(end, 2), '  Koniec', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

legend('Funkcja Celu', 'Punkt Końcowy');
grid on;
hold off;

% Wykres 3D funkcji celu
figure;
mesh(X, Y, Z);
hold on;

% Oznacz punkt końcowy na wykresie 3D
scatter3(x(end, 1), x(end, 2), f(x(end, 1), x(end, 2)), 100, 'r', 'filled');
text(x(end, 1), x(end, 2), f(x(end, 1), x(end, 2)), '  Koniec', 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');

title('3D Plot of the Objective Function');
xlabel('x1');
ylabel('x2');
zlabel('f(x1, x2)');
grid on;
hold off;

end

%% Funkcja powella

function [result, i] = Powell(f, x, d, n, Nmax, epsilon)
    syms alfa
    i = 1;
    p = [];
    while i < Nmax
        p(i, :) = x;
        for j = 1 : 1: n
            d_j = d(:, j)';
            p1 = p(i + j - 1, :) + alfa*d_j;
            alfa_j = findmin_alfa(f, p1);
            p(i+j,:) = p(i + j - 1, :) + alfa_j* d_j;
        end
        p_n = p(n + i, :);
        if abs(f(p_n(1), p_n(2)) - f(x(1), x(2))) < epsilon
            result = p(n + i, :);
            return;
        end
        for j = 1 : 1 : n-1
            d(:, j) = d(:, j+1);
        end
        d(:, n) = p(n + i, :)' - p(i, :)';
        d_n = d(:, n)';
        p_np1 = p(n + i, :) + alfa*d_n;
        alfa_n = findmin_alfa(f, p_np1);
        p(n + i + 1, :) = p(n + i, :) + alfa_n*d_n;
        x = p(n + i + 1, :);
        i = i + 1;
    end
    result = x;

end

function [alfa] = findmin_alfa(f, p1)
    syms alfa x
    g_alfa = f(p1(1), p1(2));
    fun = subs(g_alfa, alfa, x);
    fun = matlabFunction(fun);
    alfa = fminsearch(fun, 0);
end



function [val] = f_obj(param)
    
    global t;
    global G1;
    
    loc_t = t;
    
    G_test = tf([0, param(2)], [param(1), 1]);
    %G_test = G_test * pade(param(3), 3)
    
    % Calculate the step response of G1 and G_test
    step_response_G1 = step(G1, t);
    step_response_G_test = step(G_test, t);
    
    % Calculate the RMS error point-wise
    error = step_response_G1 - step_response_G_test;
    rms_error = sqrt(mean(error.^2));
    
    val = rms_error;
end
