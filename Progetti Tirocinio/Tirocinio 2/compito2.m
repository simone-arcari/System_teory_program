%% pulizia
%Per autovalori reali, distinti e negativi: gamma1=-8 e gamma2=-6
clc
clear all
close all

%% studio relativo a gamma1 e gamma2
% P_A(lambda) = det(lambda*I - A) = lambda^2 -gamma2*lambda -gamma1
% quindi gli autovalori sono:
% lambda1 = (gamma2 - sqrt(gamma2^2 + 4*gamma1))/2
% lambda1 = (gamma2 + sqrt(gamma2^2 + 4*gamma1))/2

%% caso autovalori reali e distinti
% delta = gamma2^2 + 4*gamma1 > 0
% ==>  gamma2^2 + 4*gamma1 > 0  ==>  4*gamma1 > -gamma2^2  
% ==>  gamma1 > -(gamma2^2)/4   ==>  y > -(x^2)/4
% possiamo dire che l'area sopra la parabola contiene tutti i 
% punti formati da coppie (x=gamma2, y=gamma1) tali che il delta risulta
% maggiore di zero, ovvero il caso in cui si ha autovalori reali e distinti
% si noti che i punti sulla parabola corrispondono a delta = 0 ovvero
% autovalori reali e coincidenti

figure(1)
xlim([-10, 10]);
ylim([-10, 10]);
xlabel('gamma2', 'FontSize', 16)
ylabel('gamma1', 'FontSize', 16)
title('Area blu: autovalori lamda reali e distinti', 'FontSize', 16)
hold on
grid on

x = linspace(-100, 100, 1000); % genera 1000 valori equispaziati tra -100 e 100
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

plot(x, y);

x_area = [x, fliplr(x)]; % unisce i vettori x e la versione inversa di x
y_area = [y, 100*ones(size(y))]; % unisce il vettore y e un vettore di zeri della stessa dimensione

fill(x_area, y_area, 'b');
hold off


%% caso autovalori reali e distinti entrambi negativi
% per essere negative devo coesistere le 2 seguenti condizioni:
% gamma2 - sqrt(gamma2^2 + 4*gamma1) < 0 
% gamma2 + sqrt(gamma2^2 + 4*gamma1) < 0
% la seconda condizione ingloba anche la prima quindi:
% ==>  gamma2 + sqrt(gamma2^2 + 4*gamma1) < 0  ==>  
% ==>  gamma2 < -sqrt(gamma2^2 + 4*gamma1)  ==> 
% ==>  -gamma2 > sqrt(gamma2^2 + 4*gamma1)  ==>
% ==>  (-gamma2)^2 > gamma2^2 + 4*gamma1  ==>
% ==>  gamma2^2 > gamma2^2 + 4*gamma1  ==>
% ==>  0 > 4*gamma1  ==>
% ==>  gamma1 < 0 

% non basta questo perchè essendo una disequazione irrazionale necessita
% di altre condizione:
% 1) delta > 0 ==> area sopra la parabola
% 2) gamma2 < -sqrt(gamma2^2 + 4*gamma1) < 0  == gamma2 <0

% in conclusione i punti (gamma2, gamma1) che generano autovalori reali
% e distinti negativi si trovano racchiusi tra:
% 1) gamma1 < 0
% 2) gamma2 < 0
% 3) sopra la parabola: gamma1 > -(gamma2^2)/4

figure(2)
xlim([-10, 10]);
ylim([-10, 10]);
xlabel('gamma2', 'FontSize', 16)
ylabel('gamma1', 'FontSize', 16)
title('Area verde: autovalori lamda reali e distinti entrambi negativi', 'FontSize', 16)
hold on
grid on

x = linspace(-100, 100, 1000); % genera 1000 valori equispaziati tra -100 e 100
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

plot(x, y);

x = linspace(-100, 0, 500); % genera 100 valori equispaziati tra -100 e 0
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

x_area = [x, fliplr(x)]; % unisce i vettori x e la versione inversa di x
y_area = [y, zeros(size(y))]; % unisce il vettore y e un vettore di zeri della stessa dimensione

fill(x_area, y_area, 'g');
hold off

%% caso autovalori complessi coniugati
% gli autovalori sono immagianari quando delta < 0 e ciò abbiamo già capito
% che succede per i valori (gamma2, gamma1) che stanno nell'area sottesa 
% dalla parabola

% il sottocaso di autovalori puramente immaginari si ha quando gamma2 = 0

figure(3)
xlim([-10, 10]);
ylim([-10, 10]);
xlabel('gamma2', 'FontSize', 16)
ylabel('gamma1', 'FontSize', 16)
title('Area rossa: autovalori lamda complessi coniugati', 'FontSize', 16)
hold on
grid on

x = linspace(-100, 100, 1000); % genera 1000 valori equispaziati tra -100 e 100
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

plot(x, y);

x_area = [x, fliplr(x)]; % unisce i vettori x e la versione inversa di x
y_area = [y, -100*ones(size(y))]; % unisce il vettore y e un vettore di zeri della stessa dimensione

fill(x_area, y_area, 'r');
hold off


%% caso autovalori complessi coniugati a parte reale negativa
% la parte reale è negativa se e solo se gamma2 < 0

figure(4)
xlim([-10, 10]);
ylim([-10, 10]);
xlabel('gamma2', 'FontSize', 16)
ylabel('gamma1', 'FontSize', 16)
title('Area gialla: autovalori lamda complessi coniugati a parte reale negativa', 'FontSize', 16)
hold on
grid on

x = linspace(-100, 100, 1000); % genera 100 valori equispaziati tra -100 e 100
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

plot(x, y);

x = linspace(-100, 0, 500); % genera 500 valori equispaziati tra -100 e 0
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

x_area = [x, fliplr(x)]; % unisce i vettori x e la versione inversa di x
y_area = [y, -100*ones(size(y))]; % unisce il vettore y e un vettore di zeri della stessa dimensione

fill(x_area, y_area, 'y');
hold off

%% caso autovalori reali di segno opposto
% lambda1*lambda2 < 0 (ovvero quando sono discordi)
% ==> ( (gamma2 - sqrt(gamma2^2 + 4*gamma1))/2 )*( (gamma2 + sqrt(gamma2^2 + 4*gamma1))/2 ) < 0
% ==> (gamma2^2 - gamma2^2 - 4*gamma1)/4 < 0 
% ==> (-4*gamma1)/4 < 0 
% ==> -gamma1 < 0 
% ==> gamma1 > 0 

% quando gamma1 > 0 ho i lambda discordi
% quando gamma1 < 0 ho i lambda concordi

% unendo le informazioni già ricavate si dimostra la seguente tabella che
% descrive il segno della parte reale dei lampda:
%   ___________________________________________   
%   |             | gamma_2 > 0 | gamma_2 < 0 |
%   -------------------------------------------
%   | gamma_1 < 0 |    (+|-)    |    (-|+)    |
%   -------------------------------------------
%   | gamma_1 < 0 |    (-|-)    |    (+|+)    |
%   -------------------------------------------

%  dove si intende la coppia (segno_lambda1|segno_lambda2)

figure(5)
xlim([-10, 10]);
ylim([-10, 10]);
xlabel('gamma2', 'FontSize', 16)
ylabel('gamma1', 'FontSize', 16)
title('Area MAGENTA: autovalori lamda reali e distinti di segno discorde', 'FontSize', 16)
hold on
grid on

x = linspace(-100, 100, 1000); % genera 1000 valori equispaziati tra -100 e 100
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

plot(x, y);

x_area = [x, fliplr(x)]; % unisce i vettori x e la versione inversa di x
y_area = [zeros(size(y)), 100*ones(size(y))]; % unisce il vettore y e un vettore di zeri della stessa dimensione

fill(x_area, y_area, 'm');
hold off


%% casi riassuntivi

figure(6)
xlim([-10, 10]);
ylim([-10, 10]);
xlabel('gamma2', 'FontSize', 16)
ylabel('gamma1', 'FontSize', 16)
title('Vista globale degli autovalori lambda', 'FontSize', 16)
hold on
grid on

x = linspace(-100, 100, 1000); % genera 1000 valori equispaziati tra -100 e 100
y = -(x.^2)/4; % eleva ogni elemento del vettore x al quadrato

plot(x, y);


x = linspace(-100, 0, 500); % genera 500 valori equispaziati tra -100 e 0
y = -(x.^2)/4; 
x_area = [x, fliplr(x)];
y_area = [y, zeros(size(y))]; 
fill(x_area, y_area, 'g');

x = linspace(-100, 0, 500); % genera 500 valori equispaziati tra -100 e 0
y = -(x.^2)/4;
x_area = [x, fliplr(x)];
y_area = [y, -100*ones(size(y))];
fill(x_area, y_area, 'y');

x = linspace(0, 100, 500); % genera 500 valori equispaziati tra 0 e 100
y = -(x.^2)/4;
x_area = [x, fliplr(x)];
y_area = [y, -100*ones(size(y))];
fill(x_area, y_area, 'r');

x = linspace(0, 100, 500); % genera 500 valori equispaziati tra 0 e 100
y = -(x.^2)/4;
x_area = [x, fliplr(x)];
y_area = [y, zeros(size(y))];
fill(x_area, y_area, 'w');

x = linspace(0, 100, 500); % genera 500 valori equispaziati tra 0 e 100
y = -(x.^2)/4;
x_area = [x, fliplr(x)];
y_area = [zeros(size(y)), 100*ones(size(y))];
fill(x_area, y_area, 'c');

x = linspace(-100, 0, 500); % genera 500 valori equispaziati tra -100 e 0
y = -(x.^2)/4;
x_area = [x, fliplr(x)];
y_area = [zeros(size(y)), 100*ones(size(y))];
fill(x_area, y_area, 'm');

text(-9, -2, 'reali e distinti entrambi negativi(-|-)');
text(-5, -8, 'complessi coniugati Re(-|-)');
text(1, -8, 'complessi coniugati Re(+|+)');
text(4, -2, 'reali e distinti entrambi positivi(+|+)');
text(4, 2, 'reali e distinti lambda1(-) lambda2(+)');
text(-9, 2, 'reali e distinti lambda1(+) lambda2(-)');
hold off

%% corpo principale
clc

gamma1 = input("inserisci gamma1: ");
gamma2 = input("inserisci gamma2: ");

A = [0, 1 ; gamma1, gamma2];
B = [0; 1];

lambda = eig(A, 'vector');

fprintf("Autovalori di A:\n")
z = lambda(1);
fprintf("lambda1 = %f + %fi\n", real(z), imag(z))
z = lambda(2);
fprintf("lambda2 = %f + %fi\n", real(z), imag(z))

%% secondo il Teorema di Shannon la frequenza di campionamento f_c
% deve essere almeno il doppio della frequenza massima del segnale 
% campionato f_s 
% in questo caso ciò lo si può intendere nel seguente modo:
% se vogliamo calcolare G(t) con t che che parte da zero ed 
% arriva ad un valore finale(max_time) con passo sampling_time
% il t più piccolo dopo lo zero è t = sampling_time, quindi
% l'intervallo di tempo d_tau per ricavare G(sampling_time)
% deve essere mionore o uguale di sampling_time, in modo che 
% nell'intervallo [0, d_tau] ci sia contenuato almeno 2 volte d_tau 
% ovvero ==> d_tau <= sampling_time/2

d_tau = 0.00025;
sampling_time = 0.05;
max_time = 6;
t = 0:sampling_time:max_time;


G11 = zeros(1, 1+(max_time/sampling_time));
G12 = zeros(1, 1+(max_time/sampling_time));
G21 = zeros(1, 1+(max_time/sampling_time));
G22 = zeros(1, 1+(max_time/sampling_time));

G11_analytical = 1/96 -(1/16)*exp(-4*t) + (1/12)*exp(-6*t) -(1/32)*exp(-8*t);
G12_analytical = (1/8)*exp(-4*t) -(1/4)*exp(-6*t) +(1/8)*exp(-8*t);
G21_analytical = G12_analytical;
G22_analytical = 1/12 -(1/4)*exp(-4*t) + (2/3)*exp(-6*t) -(1/2)*exp(-8*t);

for i = 1:1:(max_time/sampling_time)+1
    tau = 0:d_tau:t(i);
    G = zeros(2,2);

    for j = 1:1:(t(i)/d_tau)
        G = G + expm(A*tau(j))*B*B'*expm(A'*tau(j))*d_tau;
        G11(i) = G(1,1);
        G12(i) = G(1,2);
        G21(i) = G(2,1);
        G22(i) = G(2,2);
    end
end

figure(1)
hold on
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(1,1)', 'FontSize', 16)
title('Funzione G11(t)', 'FontSize', 16)
plot(t, G11, '-o')
plot(t, G11_analytical, '-x')
legend('G11_calcolata', 'G11_analitica')
hold off
grid on

figure(2)
hold on
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(1,2)', 'FontSize', 16)
title('Funzione G12(t)', 'FontSize', 16)
plot(t, G12, '-o')
plot(t, G12_analytical, '-x')
legend('G12 simulata', 'G12 analitica')
hold off
grid on

figure(3)
hold on
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(2,1)', 'FontSize', 16)
title('Funzione G21(t)', 'FontSize', 16)
plot(t, G21, '-o')
plot(t, G21_analytical, '-x')
legend('G21 simulata', 'G21 analitica')
hold off
grid on

figure(4)
hold on
xlabel('tempi [t]', 'FontSize', 16)
ylabel('G(2,2)', 'FontSize', 16)
title('Funzione G22(t)', 'FontSize', 16)
plot(t, G22, '-o')
plot(t, G22_analytical, '-x')
legend('G22 simulata', 'G22 analitica')
hold off
grid on