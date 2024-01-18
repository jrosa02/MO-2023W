clear all;
f = @(x) x(1)-(4.*(x(1)+x(2)).*exp(-x(1).^2 - x(2).^2)+0.1).^2
f_due = @(x1, x2) x1-(4.*(x1+x2).*exp(-x1.^2 - x2.^2)+0.1).^2
f_wrap = @(x) f_due(x(1, :), x(2, :))

[X,Y] = meshgrid(-0.6:.05:0.6);
Z = f_due(X, Y) %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! f_due
%mesh(X,Y,Z)
grid on;


P1= [-1,-1];

s = 0.5;  %  długość kroku
alpha = 0.5;  % zmniejszenie kroku po niepowodzeniu
epsilon = 1e-6;  % dokładność

contour(X,Y,Z, 100)
hold on;
x = HJ(f, P1, s, alpha, epsilon, 2)
plot(x(1), x(2), '+');
legend(["contour", "P1"], Location="northwest");

%% functions
function [new_x] = HJ_try(f, x, s, baza)
    for e = baza
        e = e';
        if f(x + s*e) < f(x)
            x = x + s.*e;
        else
            if f(x - s*e) < f(x)
                x = x - s.*e;
            end %if
        end %if
    end %while
    new_x = x;
end %HJ_try


function [x] = HJ(f, x, s, alpha, epsilon, dim)
    baza = eye(dim, dim);
    x_vec = [x];
    i = 2;
    while s>epsilon
        x_b = x;
        x = HJ_try(f, x_b, s, baza);
        if f(x) < f(x_b)
    
            while f(x) < f(x_b)
                xx_b = x_b;
                x_b = x;
                x = 2*x_b - xx_b;
                x = HJ_try(f, x, s, baza);
            end %while
            
            x = x_b;
        else %if
            s = alpha * s;
            
        end %if
    end %while
    x = x_b;
end %HJ_try
