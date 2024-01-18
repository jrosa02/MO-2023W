clear all
global G1
global t
G1=tf([0 1], [1 2 2.25 1.25]) % obiekt sterowania
Gk = tf([0 2], [1, 2, 2]) % model Gk
t = 0:0.01:20;

%Wykorzystując funkcje fminsearch, należy wykonać optymalizację nastaw regulatora PID tak, aby odpowiedź układu regulacji na zmianę wartości referencyjnej - skok jednostkowy była zbliżona do dynamiki modelu Gk(s) (dominująca stała czasowa τ≈1.25[s])
x0 = [0,0,0]
fun = @ident;
%------
A = []; % No other constraints
b = [];
Aeq = [];
beq = [];
lb = [0, 0, 0];
ub = [100, 100, 100];
nonlcon = [];
intcon = [];
options = optimoptions("ga",'PlotFcn', @gaplotchange, 'Display','iter');
options.InitialPopulationMatrix = x0;
[x, fval, exitflag] = ga(fun, 3, A, b, Aeq, beq, lb, ub, nonlcon, intcon, options);

reg=pid(x(1), x(2), x(3));

G1c = feedback(reg*G1, 1);
Gk = tf([0 2], [1, 2, 2]);

figure;
step(G1 ,G1c, Gk, t);
legend(["G1", "G1c", "Gk"]);
grid on;


%% funkcje

function blad = ident(x)
global G1 
global t

Gk = tf([0 2], [1, 2, 2]);

reg = pid(x(1), x(2), x(3));

G1c = feedback(reg*G1, 1);

[y1] = step(Gk, t); % wzor
[y2] = step(G1c, t); % dopasowanie

e = y1-y2;
blad = sum(e.^2);

end %ident

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

