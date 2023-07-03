%%Scegliere un tempo t_bar (10 sec) per ognuno dei 4 casi gli
%%Diamo un ingresso e lo simuliamo con Lsim (non facciamo nessun calcolo)
%%Ricordarsi di capire come funziona Lsim
%%Ricordarsi di salvare il valore di t_bar scelto
%%Salvarsi i grafici prodotti dopo questa cosa
%%Gli ingressi dobbiamo sceglierli in maniera variegata in base anche al
%%caso in esame (stabili/instabili)

t_bar = 3;
t = 0:sampling_time:t_bar;

input("\nPremi un tasto per andare avanti col punto2:")

C = eye(2);
D = 0;

x0 = [0; 0];

sys = ss(A,B,C,D);

%% gamma1: -8, gamma2: -6, --> lambda1: -2, lambda2: -4
if gamma1 == -8 && gamma2 == -6
    u_0 = sin(t);
    state = lsim(sys, u_0, t, x0);
    x_bar = state(end, :);
end

%% gamma1: 2, gamma2: -1, --> lambda1: 1, lambda2: -2
if gamma1 == 2 && gamma2 == -1
    u_0 = sin(t);
    state = lsim(sys, u_0, t, x0);
    x_bar = state(end, :);
end

%% gamma1: -2, gamma2: -2, --> lambda1: -1+j, lambda2: -1-j
if gamma1 == -2 && gamma2 == -2
    u_0 = sin(t);
    state = lsim(sys, u_0, t, x0);
    x_bar = state(end, :);
end

%% gamma1: -2, gamma2: 2, --> lambda1: 1-j, lambda2: 1+j
if gamma1 == -2 && gamma2 == 2
    u_0 = sin(t);
    state = lsim(sys, u_0, t, x0);
    x_bar = state(end, :);
end

%% Grafici
figure(5)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(t, state(:,1), '-x')
hold on
plot(t, state(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

punto3