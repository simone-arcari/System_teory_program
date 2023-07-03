% Insieme di tempi T1, T2, ..., Tm (incluso t_bar)
% Per ognuno di questi calcolare un controllo u, partendo da x0 fino a x_bar

input("\nPremi un tasto per andare avanti col punto3:")

T_i = 2:1:10;

x(1,1) = 1%x_bar(1,1);
x(2,1) = 2%x_bar(1,2);

u1 = zeros(1, 1 + T_i(1)/sampling_time);
u2 = zeros(1, 1 + T_i(2)/sampling_time);
u3 = zeros(1, 1 + T_i(3)/sampling_time);
u4 = zeros(1, 1 + T_i(4)/sampling_time);
u5 = zeros(1, 1 + T_i(5)/sampling_time);
u6 = zeros(1, 1 + T_i(6)/sampling_time);
u7 = zeros(1, 1 + T_i(7)/sampling_time);
u8 = zeros(1, 1 + T_i(8)/sampling_time);
u9 = zeros(1, 1 + T_i(9)/sampling_time);

tau1 = 0:sampling_time:T_i(1);
tau2 = 0:sampling_time:T_i(2);
tau3 = 0:sampling_time:T_i(3);
tau4 = 0:sampling_time:T_i(4);
tau5 = 0:sampling_time:T_i(5);
tau6 = 0:sampling_time:T_i(6);
tau7 = 0:sampling_time:T_i(7);
tau8 = 0:sampling_time:T_i(8);
tau9 = 0:sampling_time:T_i(9);



for i = 1:1:1%length(T_i)

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


    beta = pinv(G)*x;

    tau = 0:sampling_time:T_i(i);
    
    for j = 1:1:length(tau)
        switch i
            case 1
                u1(j) = B'*exp(A'*(T_i(1)-tau(j)))*beta;

            case 2
                u2(j) = B'*exp(A'*(T_i(2)-tau(j)))*beta;

            case 3
                u3(j) = B'*exp(A'*(T_i(3)-tau(j)))*beta;

            case 4
                u4(j) = B'*exp(A'*(T_i(4)-tau(j)))*beta;

            case 5
                u5(j) = B'*exp(A'*(T_i(5)-tau(j)))*beta;

            case 6
                u6(j) = B'*exp(A'*(T_i(6)-tau(j)))*beta;

            case 7
                u7(j) = B'*exp(A'*(T_i(7)-tau(j)))*beta;

            case 8
                u8(j) = B'*exp(A'*(T_i(8)-tau(j)))*beta;

            case 9
                u9(j) = B'*exp(A'*(T_i(9)-tau(j)))*beta;

            otherwise
                disp('errore');
        end  
    end

  
end

%% Grafici Ingressi u_i
close all

figure(1)
plot(tau1, u1)
grid on

figure(2)
plot(tau2, u2)
grid on

figure(3)
plot(tau3, u3)
grid on

figure(4)
plot(tau4, u4)
grid on

figure(5)
plot(tau5, u5)
grid on

figure(6)
plot(tau6, u6)
grid on

figure(7)
plot(tau7, u7)
grid on

figure(8)
plot(tau8, u8)
grid on

figure(9)
plot(tau9, u9)
grid on

%% 
stato1 = lsim(sys, u1, tau1, x0);
stato2 = lsim(sys, u2, tau2, x0);
stato3 = lsim(sys, u3, tau3, x0);
stato4 = lsim(sys, u4, tau4, x0);
stato5 = lsim(sys, u5, tau5, x0);
stato6 = lsim(sys, u6, tau6, x0);
stato7 = lsim(sys, u7, tau7, x0);
stato8 = lsim(sys, u8, tau8, x0);
stato9 = lsim(sys, u9, tau9, x0);


figure(10)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau1, stato1(:,1), '-x')
hold on
plot(tau1, stato1(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(11)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau2, stato2(:,1), '-x')
hold on
plot(tau2, stato2(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(12)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau3, stato3(:,1), '-x')
hold on
plot(tau3, stato3(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(13)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau4, stato4(:,1), '-x')
hold on
plot(tau4, stato4(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(14)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau5, stato5(:,1), '-x')
hold on
plot(tau5, stato5(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(15)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau6, stato6(:,1), '-x')
hold on
plot(tau6, stato6(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(16)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau7, stato7(:,1), '-x')
hold on
plot(tau7, stato7(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(17)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau8, stato8(:,1), '-x')
hold on
plot(tau8, stato8(:,2), '-x')
legend('x1', 'x2')
hold off
grid on

figure(18)
xlabel('tempi [t]', 'FontSize', 16)
ylabel('X', 'FontSize', 16)
title('Andamento Stati', 'FontSize', 16)
plot(tau9, stato9(:,1), '-x')
hold on
plot(tau9, stato9(:,2), '-x')
legend('x1', 'x2')
hold off
grid on









G

Ginversa = pinv(G)

beta = pinv(G)*x


