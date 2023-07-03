% Insieme di tempi T1, T2, ..., Tm (incluso t_bar)
% Per ognuno di questi calcolare un controllo u, partendo da x0 fino a x_bar

input("\nPremi un tasto per andare avanti col punto3:")

T_i = 2:1:10;

x_f(1,1) = x_bar(1,1);
x_f(2,1) = x_bar(1,2);

u1 = [];
u2 = [];
u3 = [];
u4 = [];
u5 = [];
u6 = [];
u7 = [];
u8 = [];
u9 = [];

clear tau
syms tau

tau1 = 0:sampling_time:T_i(1);
tau2 = 0:sampling_time:T_i(2);
tau3 = 0:sampling_time:T_i(3);
tau4 = 0:sampling_time:T_i(4);
tau5 = 0:sampling_time:T_i(5);
tau6 = 0:sampling_time:T_i(6);
tau7 = 0:sampling_time:T_i(7);
tau8 = 0:sampling_time:T_i(8);
tau9 = 0:sampling_time:T_i(9);



for i = 1:1:length(T_i)

    if gamma1 == -8 && gamma2 == -6
        G(1,1) = 1/96 + (-1/16)*exp(-4*T_i(i)) + (1/12)*exp(-6*T_i(i)) + (-1/32)*exp(-8*T_i(i));
        G(1,2) = (1/8)*exp(-4*T_i(i)) + (-1/4)*exp(-6*T_i(i)) + (1/8)*exp(-8*T_i(i));
        G(2,1) = G(1,2);
        G(2,2) = 1/12 + (-1/4)*exp(-4*T_i(i)) + (2/3)*exp(-6*T_i(i)) + (-1/2)*exp(-8*T_i(i));
    end

    if gamma1 == 2 && gamma2 == -1
        G(1,1) = (-1/36)*exp(-4*T_i(i)) + (2/9)*exp(-T_i(i)) + (1/18)*exp(2*T_i(i)) -25/36;
        G(1,2) = (1/18)*exp(-4*T_i(i)) + (-1/9)*exp(-T_i(i)) + (1/18)*exp(2*T_i(i));
        G(2,1) = G(1,2);
        G(2,2) = (-1/9)*exp(-4*T_i(i)) + (-4/9)*exp(-T_i(i)) + (1/18)*exp(2*T_i(i)) +1/2;
    end
    
    if gamma1 == -2 && gamma2 == -2
        G(1,1) = (1/8)*(exp(-2*T_i(i))).*(-sin(2*T_i(i))+cos(2*T_i(i))) + (-1/4)*exp(-2*T_i(i)) + (1/8);
        G(1,2) = (-1/4)*(exp(-2*T_i(i))).*cos(2*T_i(i)) + (1/4)*exp(-2*T_i(i));
        G(2,1) = G(1,2);
        G(2,2) = (1/4)*(exp(-2*T_i(i))).*(sin(2*T_i(i))+cos(2*T_i(i))) + (-1/2)*exp(-2*T_i(i)) + (1/4);
    end
    
    if gamma1 == -2 && gamma2 == 2
        G(1,1) = (-1/8)*(exp(2*T_i(i))).*(sin(2*T_i(i))+cos(2*T_i(i))) + (1/4)*exp(2*T_i(i)) + (-1/8);
        G(1,2) = (-1/4)*(exp(2*T_i(i))).*cos(2*T_i(i)) + (1/4)*exp(2*T_i(i));
        G(2,1) = G(1,2);
        G(2,2) = (1/4)*(exp(2*T_i(i))).*(sin(2*T_i(i))-cos(2*T_i(i))) + (1/2)*exp(2*T_i(i)) + (-1/4);
    end


    beta = pinv(G)*x_f;
    u_cappuccio = B'*expm(A'*(T_i(i)-tau))*beta;

    time = 0:sampling_time:T_i(i);
    U = [];
    
    for k = 1:1:length(time)
        switch i
            case 1
                u1 = [u1 double(subs(u_cappuccio, tau, time(k)))];

            case 2
                u2 = [u2 double(subs(u_cappuccio, tau, time(k)))];

            case 3
                u3 = [u3 double(subs(u_cappuccio, tau, time(k)))];

            case 4
                u4 = [u4 double(subs(u_cappuccio, tau, time(k)))];

            case 5
                u5 = [u5 double(subs(u_cappuccio, tau, time(k)))];

            case 6
                u6 = [u6 double(subs(u_cappuccio, tau, time(k)))];

            case 7
                u7 = [u7 double(subs(u_cappuccio, tau, time(k)))];

            case 8
                u8 = [u8 double(subs(u_cappuccio, tau, time(k)))];

            case 9
                u9 = [u9 double(subs(u_cappuccio, tau, time(k)))];

            otherwise
                disp('errore');
        end   
   end

  
end

stato1 = lsim(sys, u1, tau1, x0);
stato2 = lsim(sys, u2, tau2, x0);
stato3 = lsim(sys, u3, tau3, x0);
stato4 = lsim(sys, u4, tau4, x0);
stato5 = lsim(sys, u5, tau5, x0);
stato6 = lsim(sys, u6, tau6, x0);
stato7 = lsim(sys, u7, tau7, x0);
stato8 = lsim(sys, u8, tau8, x0);
stato9 = lsim(sys, u9, tau9, x0);

%% Grafici 

figure(3)
hold on
xlabel('tempi [t]', 'FontSize', 16)
ylabel('U', 'FontSize', 16)
title('Andamento Ingressi per vari T_i', 'FontSize', 16)
plot(tau1, u1)
plot(tau2, u2)
plot(tau3, u3)
plot(tau4, u4)
plot(tau5, u5)
plot(tau6, u6)
plot(tau7, u7)
plot(tau8, u8)
plot(tau9, u9)
grid on
 

figure(4)
hold on
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)

plot(tau1, stato1(:,1), 'color', 'red')
plot(tau1, stato1(:,2), 'color', 'blue')

plot(tau2, stato2(:,1), 'color', 'red')
plot(tau2, stato2(:,2), 'color', 'blue')

plot(tau3, stato3(:,1), 'color', 'red')
plot(tau3, stato3(:,2), 'color', 'blue')

plot(tau4, stato4(:,1), 'color', 'red')
plot(tau4, stato4(:,2), 'color', 'blue')

plot(tau5, stato5(:,1), 'color', 'red')
plot(tau5, stato5(:,2), 'color', 'blue')

plot(tau6, stato6(:,1), 'color', 'red')
plot(tau6, stato6(:,2), 'color', 'blue')

plot(tau7, stato7(:,1), 'color', 'red')
plot(tau7, stato7(:,2), 'color', 'blue')

plot(tau8, stato8(:,1), 'color', 'red')
plot(tau8, stato8(:,2), 'color', 'blue')

plot(tau9, stato9(:,1), 'color', 'red')
plot(tau9, stato9(:,2), 'color', 'blue')

grid on




