%maksymalna objętość kuli przy minimalnej powierzchni stożka

clear;
ball = @(R) 4/3*pi*(R^2)
spike = @(r, b) pi*r*(r + b)
p_eval = @(a, b) (a+2*b)/2

fun = @(x)((spike(x(1), x(2))/ball(R_eval(x(1), x(2)))))
A = [1, -1];
b = 0;

x0 = [20,90]; %[r, l]
lb = [0.001, 0.001];
ub = [100,100];

options = optimoptions("ga",'PlotFcn', @gaplotchange, 'Display','iter');
options.InitialPopulationMatrix = x0;
[x,fval,exitflag,output,population,scores] = ga(fun, 2, A, b, [], [], lb, ub, [], [], options)

x(2)/x(1)
ball(R_eval(x(1), x(2)))

function [R] = R_eval(a, b)
    p_eval = @(a, b)[(a+2*b)/2];
    p = p_eval(a, b);
    R = sqrt(((p-a)*(p-b)^2)/p);
end %function

function state = gaplotchange(options, state, flag)
persistent last_best
if(strcmp(flag,'init'))
xlim([1,options.MaxGenerations]);
axx = gca;
axx.YScale = 'log';
hold on;
xlabel Generation
X = sprintf('Best value: %.5f', min(state.Score));
title(X)
end
X = sprintf('Best value: %.5f', min(state.Score));
title(X)
best = min(state.Score);
plot(state.Generation,best,'o', 'MarkerFaceColor',[1 0.5 0], ...
'MarkerEdgeColor', [1 0.5 0]);
end
