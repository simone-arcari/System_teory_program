%% Inizio del programma
clc
clear all


fprintf("**************************READ ME**************************\n\n")

fprintf("\tCaso 1) Autovalori Reali Negativi\n")
fprintf("\t\t\tgamma1: -8\n")
fprintf("\t\t\tgamma2: -6\n\n")
fprintf("\t\t\tlambda1: -2\n")
fprintf("\t\t\tlambda2: -4\n\n")

fprintf("\tCaso 2) Autovalori Reali Uno Positivo e Uno Negativi\n")
fprintf("\t\t\tgamma1: 2\n")
fprintf("\t\t\tgamma2: -1\n\n")
fprintf("\t\t\tlambda1: 1\n")
fprintf("\t\t\tlambda2: -2\n\n")

fprintf("\tCaso 3) Autovalori Complessi Coniugati Stabili\n")
fprintf("\t\t\tgamma1: -2\n")
fprintf("\t\t\tgamma2: -2\n\n")
fprintf("\t\t\tlambda1: -1+j\n")
fprintf("\t\t\tlambda2: -1-j\n\n")

fprintf("\tCaso 4) Autovalori Complessi Coniugati Instabili\n")
fprintf("\t\t\tgamma1: -2\n")
fprintf("\t\t\tgamma2: 2\n\n")
fprintf("\t\t\tlambda1: 1+j\n")
fprintf("\t\t\tlambda2: 1-j\n\n")

fprintf("***********************************************************\n\n")

gamma1 = input("inserisci gamma1: ");
gamma2 = input("inserisci gamma2: ");

A = [0, 1 ; gamma1, gamma2];
B = [0; 1];

lambda = eig(A, 'vector');

fprintf("\nAutovalori di A:\n")
z = lambda(1);
fprintf("lambda1 = %f + %fi\n", real(z), imag(z))
z = lambda(2);
fprintf("lambda2 = %f + %fi\n", real(z), imag(z))

%% In questa sezione si calcolano le matrici G(0), G(sampling_time), 
% G(2*sampling_time), G(3*sampling_time), ... , G(max_time)
%
% Per fare ciò simuliamo il calcolo integrale con la definizione di
% integrale di Riemann:
%
% ESEMPIO 1: G(t = sampling_time) = Funzione_integranda(0)*d_tau + 
% Funzione_integranda(d_tau)*d_tau + Funzione_integranda(2*d_tau)*d_tau + 
% + ... + Funzione_integranda(sampling_time)*d_tau 
%
% ESEMPIO 2: G(t = 3*sampling_time) = Funzione_integranda(0)*d_tau + 
% Funzione_integranda(d_tau)*d_tau + Funzione_integranda(2*d_tau)*d_tau + 
% + ... + Funzione_integranda(3*sampling_time)*d_tau
% 
% Secondo il Teorema di Shannon la frequenza di campionamento f_c
% deve essere almeno il doppio della frequenza massima del segnale 
% campionato f_s
% In questo caso ciò lo si può intendere nel seguente modo:
% se vogliamo calcolare G(t) con t che che parte da zero ed 
% arriva ad un valore finale(max_time) con passo sampling_time
% il t più piccolo dopo lo zero è t = sampling_time, quindi
% l'intervallo di tempo d_tau per ricavare G(t=sampling_time)
% deve essere mionore o uguale di sampling_time, in modo che 
% nell'intervallo [0, sampling_time] ci sia contenuto almeno 2 volte d_tau 
% ovvero ==> d_tau <= sampling_time/2
%
% Facendo tendere a zero d_tau e sampling_time si otterà un integrazione
% sempre più precisa ma con un aumento notevole del costo computazionale

d_tau = 0.00025;
sampling_time = 0.025;
max_time = 20;
t = 0:sampling_time:max_time;


G11 = zeros(1, 1+(max_time/sampling_time));
G12 = zeros(1, 1+(max_time/sampling_time));
G21 = zeros(1, 1+(max_time/sampling_time));
G22 = zeros(1, 1+(max_time/sampling_time));


contatore = 1;
massimo_contatore = (max_time/sampling_time);
fprintf('\n-->Percentuale completamento del calcolo integrale: 0######')

G = zeros(2,2);
for i = 2:1:(max_time/sampling_time)+1
    tau = t(i-1):d_tau:t(i);
    %M = (B')*expm(A'*tau(i));


    % codice per stampare la parcentuale di completamento del calcolo------
    percentuale = 100*(contatore/massimo_contatore);

    if percentuale < 10
        fprintf('\b\b\b\b\b\b')
    else
        fprintf('\b\b\b\b\b\b\b')
    end
    
    fprintf('%.2f%%\n', 100*(contatore/massimo_contatore))
    contatore = contatore + 1;
    %----------------------------------------------------------------------

    % Sommatoria di Riemann
    for j = 1:1:(sampling_time/d_tau)
        M1 = expm(A*tau(j))*B;
        M2 = (B')*expm(A'*tau(j));

        G = G + M1*M2*d_tau;
        %%beta = inv(M2)*u()

        G11(i) = G(1,1);
        G12(i) = G(1,2);
        G21(i) = G(2,1);
        G22(i) = G(2,2);
    end
end

%% gamma1: -8, gamma2: -6, --> lambda1: -2, lambda2: -4
if gamma1 == -8 && gamma2 == -6
    G11_analytical = 1/96 + (-1/16)*exp(-4*t) + (1/12)*exp(-6*t) + (-1/32)*exp(-8*t);
    G12_analytical = (1/8)*exp(-4*t) + (-1/4)*exp(-6*t) + (1/8)*exp(-8*t);
    G21_analytical = G12_analytical;
    G22_analytical = 1/12 + (-1/4)*exp(-4*t) + (2/3)*exp(-6*t) + (-1/2)*exp(-8*t);
end

%% gamma1: 2, gamma2: -1, --> lambda1: 1, lambda2: -2
if gamma1 == 2 && gamma2 == -1
    G11_analytical = (-1/36)*exp(-4*t) + (2/9)*exp(-t) + (1/18)*exp(2*t) -25/36;
    G12_analytical = (1/18)*exp(-4*t) + (-1/9)*exp(-t) + (1/18)*exp(2*t);
    G21_analytical = G12_analytical;
    G22_analytical = (-1/9)*exp(-4*t) + (-4/9)*exp(-t) + (1/18)*exp(2*t) +1/2;
end

%% gamma1: -2, gamma2: -2, --> lambda1: -1+j, lambda2: -1-j
if gamma1 == -2 && gamma2 == -2
    G11_analytical = (1/8)*(exp(-2*t)).*(-sin(2*t)+cos(2*t)) + (-1/4)*exp(-2*t) + (1/8);
    G12_analytical = (-1/4)*(exp(-2*t)).*cos(2*t) + (1/4)*exp(-2*t);
    G21_analytical = G12_analytical;
    G22_analytical = (1/4)*(exp(-2*t)).*(sin(2*t)+cos(2*t)) + (-1/2)*exp(-2*t) + (1/4);
end

%% gamma1: -2, gamma2: 2, --> lambda1: 1-j, lambda2: 1+j
if gamma1 == -2 && gamma2 == 2
    G11_analytical = (-1/8)*(exp(2*t)).*(sin(2*t)+cos(2*t)) + (1/4)*exp(2*t) + (-1/8);
    G12_analytical = (-1/4)*(exp(2*t)).*cos(2*t) + (1/4)*exp(2*t);
    G21_analytical = G12_analytical;
    G22_analytical = (1/4)*(exp(2*t)).*(sin(2*t)-cos(2*t)) + (1/2)*exp(2*t) + (-1/4);
end

%% Grafici
close all

figure(1)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(1,1)', 'FontSize', 16)
title('Funzione G11(t)', 'FontSize', 16)
plot(t, G11, '-o')
hold on
plot(t, G11_analytical, '-x')
legend('G11 calcolata', 'G11 analitica')
hold off
grid on

figure(2)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(1,2)', 'FontSize', 16)
title('Funzione G12(t)', 'FontSize', 16)
plot(t, G12, '-o')
hold on
plot(t, G12_analytical, '-x')
legend('G12 simulata', 'G12 analitica')
hold off
grid on

figure(3)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(2,1)', 'FontSize', 16)
title('Funzione G21(t)', 'FontSize', 16)
plot(t, G21, '-o')
hold on
plot(t, G21_analytical, '-x')
legend('G21 simulata', 'G21 analitica')
hold off
grid on

figure(4)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(2,2)', 'FontSize', 16)
title('Funzione G22(t)', 'FontSize', 16)
plot(t, G22, '-o')
hold on
plot(t, G22_analytical, '-x')
legend('G22 simulata', 'G22 analitica')
hold off
grid on

punto2