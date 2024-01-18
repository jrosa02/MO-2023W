clear all;
syms x;
syms y;

%tutaj funkcja
f = x .* exp(-x.^2 - y.^2);
%f= 0.3*x+0.1*y+(-3.5+0.5*x.^2+0.5*y.^2).^2+100*x .* exp(-x.^2 - y.^2);

grad = gradient(f, [x, y])

sol = solve(grad, [x, y])

hes = hessian(f, [x, y])

for i = 1:length(sol.x)
    double([sol.x(i), sol.y(i)])
    solll = subs(hes, [x, y], [sol.x(i), sol.y(i)]);
    all(0<double(eig(solll))) % if positive definite
    all(0>double(eig(solll))) % if negative definite
end
fsurf(f)