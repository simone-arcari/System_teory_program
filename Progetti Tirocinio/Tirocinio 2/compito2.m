clc
clear all
close all

gamma1 = -1;
gamma2 = 2;

A = [0, 1 ; gamma1, gamma2];
B = [0; 1];

sampling_time = 1;
max_time = 3;
t = 0:sampling_time:max_time;
d_tau = 0.1;

G11 = zeros(1,max_time+1);

for i = 1:1:max_time+1
    tau = 0:d_tau:t(i);
    G = zeros(2,2);

    for j = 1:1:t(i)+1
        G = G + expm(-A*tau(j))*B*B'*expm(-A'*tau(j));
        G11(i) = G(1,1);
    end
end

figure(1)
plot(t, G11)