syms x1;
syms x2;
f_syms = 0.3*x1+0.1*x2+(-3.5+0.5*x1.^2+0.5*x2.^2).^2+100*x1 .* exp(-x1.^2 - x2.^2)
f= @(x1, x2) 0.3*x1+0.1*x2+(-3.5+0.5*x1.^2+0.5*x2.^2).^2+100*x1 .* exp(-x1.^2 - x2.^2)

figure;
[X, Y] = meshgrid(-3:0.1:3, -3:0.1:3);
Z = f(X, Y);
mesh(X, Y, Z);
contour(X, Y, Z, 40);

n = 5;
peak = [ 0.8, 0];
pitt = [-0.8, 0];

o = f(0.8, 0)
o = f(-0.8, 0)

f = @(points) distance(points, peak, pitt);

x0 =  zeros(n, 2);
lin = linspace(peak(1), pitt(1), n+2);
x0(:, 1) = lin(:, 2:end-1);
lin = linspace(peak(2), pitt(2), n+2);
x0(:, 2) = lin(:, 2:end-1);
x0 = reshape(x0, [n * 2, 1]);

orginal = f(x0)

x = fminsearch(f, x0);

fmin_search = f(x)
x = reshape(x, [n,2])

options = optimset('Display', 'iter');

options.InitialPopulationMatrix = x0;
[x, ~] = ga(f, n*2, [], [], [], [], [], [], [], options);
ga_search = f(x)
x = reshape(x, [n,2]) %distanceeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee

f = @(x1, x2) 0.3*x1 + 0.1*x2 + (-3.5 + 0.5*x1.^2 + 0.5*x2.^2).^2 + 100*x1 .* exp(-x1.^2 - x2.^2);
[X, Y] = meshgrid(-3:0.1:3, -3:0.1:3);
Z = f(X, Y);
figure
mesh(X, Y, Z);

clf
figure;
f = @(x1, x2) 0.3*x1 + 0.1*x2 + (-3.5 + 0.5*x1.^2 + 0.5*x2.^2).^2 + 100*x1 .* exp(-x1.^2 - x2.^2);
[X, Y] = meshgrid(-3:0.1:3, -3:0.1:3);
Z = f(X, Y);
hold on
mesh(X, Y, Z);
for i = 2:n
    plot3(x(i, 1), x(i, 2), f(x(i, 1), x(i, 2)),'-o');
end
hold off



function [dis] = distance(points_ar, peak, pitt)
    f = @(x1, x2) 0.3*x1 + 0.1*x2 + (-3.5 + 0.5*x1.^2 + 0.5*x2.^2).^2 + 100*x1 .* exp(-x1.^2 - x2.^2);
    F = @(x) f(x(1), x(2));
    n = 5;

    points = reshape(points_ar, [n,2]);


    % Calculate the distance without penalty term
    dis = sqrt((peak(1) - points(1, 1))^2 + (peak(2) - points(1, 2))^2 + (F(peak) - F(points(1, :)))^2);
    desired_spacing = dis;
    spacing_penalty = 0;

    for i = 2:n
        last_dist = sqrt((points(i, 1) - points(i-1, 1))^2 + (points(i, 2) - points(i-1, 2))^2 + (F(points(i, :)) - F(points(i-1, :)))^2);
        dis = dis + last_dist;
        spacing_penalty = spacing_penalty + abs(last_dist - desired_spacing);
    end
    
    last_dist = sqrt((pitt(1) - points(n, 1))^2 + (pitt(2) - points(n, 2))^2 + (F(pitt) - F(points(n, :)))^2);
    dis = dis + last_dist;
    spacing_penalty = spacing_penalty + abs(last_dist - desired_spacing);
    
    
    % Adjust the weight of the penalty term based on your preference
    penalty_weight = 0; % Adjust this value based on your preference

    dis = dis + penalty_weight * spacing_penalty;
end

