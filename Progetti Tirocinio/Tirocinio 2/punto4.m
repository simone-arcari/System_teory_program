input("\nPremi un tasto per andare avanti col punto4:")

J_t_bar = int(1000*sin(tau)*1000*sin(tau), tau, 0, t_bar);

%% Grafici
figure(6)

hold on
xlabel('T_i', 'FontSize', 16)
ylabel('J(T_i)', 'FontSize', 16)
title('Costi controlli ottimi', 'FontSize', 16)
plot(T_i, J_Ti, '-o')
grid on

figure(7)
xlabel('t_bar', 'FontSize', 16)
ylabel('J(t_bar)', 'FontSize', 16)
title('Costo controllo u_0', 'FontSize', 16)
plot(t_bar, J_t_bar, '-o')
grid on