%Zadanie 2) Optymalizacja regulatora PID
%Na podstawie pomiarów odpowiedzi układu sterowania na skok jednostkowy wartości referencyjnej należy wykonać optymalizację nastaw regulatora wg krzywej narastania Gk(s) wartości wyjścia z obiektu sterowanego, wykorzystując optymalizację:
%Warunek początkowy x0 = [0,0,0],  gdzie x0 = [Kp,Ki,Kd], estymowane Kp=x(1), Ki=x(2),  Kd=x(3)
%Należy wykorzystać funkcję Matlaba C = pid(Kp,Ki,Kd) równanie reg w formie operatorowej
%gdzie:
clear all
global G1
global t
G1=tf([0 1], [1 2 2.25 1.25]) % obiekt sterowania
Gk = tf([0 2], [1, 2, 2]) % model Gk
t = 0:0.01:20;
%C = pid(x(1),x(2),x(3))         % optymalizowany kontroler
%a)    Optymalizacja bez ograniczeń
%Wykorzystując funkcje fminsearch, należy wykonać optymalizację nastaw regulatora PID tak, aby odpowiedź układu regulacji na zmianę wartości referencyjnej - skok jednostkowy była zbliżona do dynamiki modelu Gk(s) (dominująca stała czasowa τ≈1.25[s])
%                       2
%Gk(s) = -------------------
%              s^2 + 2 s + 2
options = optimset();
x0 = [1, 1, 1];
[x,fval, exitflag] = fminsearch(@ident, x0, options)

reg=pid(x(1), x(2), x(3));

G1c = feedback(reg*G1, 1);
Gk = tf([0 2], [1, 2, 2]);

figure;
step(G1 ,G1c, Gk, t);
legend(["G1", "G1c", "Gk"]);
grid on;

%% b)    Optymalizacja z ograniczeniami
%Wykorzystując funkcje fmincon, należy wykonać optymalizację nastaw regulatora PID tak aby odpowiedz układu regulacji na zmianę wartości referencyjnej - skok jednostkowy była zbliżona do dynamiki modelu Gk(s) (dominująca stała czasowa τ≈1.25[s])
%przy ograniczeniach na nastawy P,I,D, regulator rzeczywisty i P∈[0.2…5], I∈[0.2…10], D∈[0.2…10]
%przy ograniczeniach równościowych i  nierównościowych spełniających kryterium Hurwitza dla transmitancji zastępczej układu regulator-obiekt.
%Wskazówka:
%Rozważyć rożne warunki początkowe np. x0 = [1,1,1],.
x0 = [1, 1, 1];
opts = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'sqp');
x = fmincon(@ident ,x0, [], [], [], [], [], [], @ogr,opts)

reg=pid(x(1), x(2), x(3));

G1c = feedback(reg*G1, 1);
Gk = tf([0 2], [1, 2, 2]);

figure;
step(G1, G1c, Gk, t);
legend(["G1", "G1c", "Gk"]);
grid on;

%% functions

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

function [c, ceq] = Hurwitz(x)
    global G1 
    c = [0, 0];
    ceq = [];
    
    reg = pid(x(1), x(2), x(3));
    G1c = feedback(reg*G1, 1);
    
    a_vec = G1c.Denominator{1};
    
    if all(a_vec > 0)
        %do nothing
    elseif all(a_vec < 0)
        a_vec = -a_vec;
    else
        %does not satisify 1st condition
        c = [abs(min(a_vec)), abs(min(a_vec))]; 
        return;
    end %if
    
    [H, delta] = hurwitz_matrix(a_vec);
    
    if all(delta > 0)
        c = [-1, -1];
    else
        c = [1, 1];
    end %if

end %Hurwitz

function [c, ceq] = ogr(x)
    ceq = [];
    if 0.2<x(1) && x(1)<5 && 0.2<x(2) && x(2)<10 && 0.2<x(3) && x(3)<10
        [c, ceq] = Hurwitz(x);
    else
        c = [1, 1];
    end %if
end %ogr

%% hurwitz stability where delta is determinnats ofminor vector

function [H, delta] = hurwitz_matrix(p)

    n = numel(p) - 1;
    p1 = p(2:2:end);
    p2 = p(1:2:end);
    
    if isnumeric(p)
        H = zeros(n, n);
        delta = zeros(n, 1);
    else
        H = zeros(n, n);
        delta = sym(zeros(n, 1));
    end
    
    i = 0;
    for k = 1:n
        if mod(k, 2)
            H(k, i+[1:numel(p1)]) = p1;
        else
            H(k, i+[1:numel(p2)]) = p2;
        end
    end
    
    for k = 1:n
        delta(k) = det(H(1:k, 1:k));
    end

end
