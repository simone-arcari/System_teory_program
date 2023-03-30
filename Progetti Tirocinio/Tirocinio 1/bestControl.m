function [u, P, w] = bestControl(A,B, nu, x_bar)
%   BESTCONTROL: This function calculates the optimal control u, assuming x(0)=0
%   Step 1: Calculate P_nu
%   Step 2: Calculate G = P*P' , find beta so that G*beta = x_bar and define w = P'*beta
%   Step 3/4: Break up w into nu subvectors and assigned them as u(h) = u_h for h=0,..,nu-1

%STEP1 ********************************************************************
    
    P = B;                                          %Setting up P_nu as B
    
    for i = 1:1:nu-1                                %Index starts at 1 because P_nu is already setted up as B

        M = (A^i)*B;                                %Calculate A^i*B
        P = [P M];                                  %Every iteration adds matrix A^i*B to current P_nu
    end

%STEP2 ********************************************************************

    P_t = P';                                       %Calculate P traspose 
    G = P*P_t;                                      %Calculate G
   
    G_inv = pinv(G);                                %Invert G 
    beta = G_inv*x_bar;                             %Calculate beta

    x = G*beta;                                     %vector x used to check beta validity
   
    epsilon = ones(size(x_bar, 1), size(x_bar, 2) ) * 0.0001;     %Thresold vector to confront x_bar and x made 
        

    if abs(x-x_bar) > epsilon                       %if G*beta is diffrent from x_bar beta is an incorrect solution

        fprintf("[MESSAGE] --> La soluzione beta trovata per k = %d non è ammissibile\n", nu)
        u = -1;
        P = -1;
        w = -1;
        return
    end

    w = P_t*beta;                                   %Calculate w

%STEP3/4 ******************************************************************
    
    p = size(B, 2);                                 %Return number of coloums of B

    topCut = 1 + p*(nu-1);                          %Calculate top index to partition w
    bottomCut = topCut + (p-1);                     %Calculate bottom index to partition w
    u = w(topCut : bottomCut);                      %First w partition 
    
    for i = 1:1:nu-1                                %For loop to partition w; Starting from 1 because we already calculate u_0
    
        topCut = topCut - p;                        %Top index progression 
        bottomCut = topCut + (p -1);                %Bottom index progression
        u = [u; w(topCut : bottomCut)];             %Every iteration adds u_i to the already existing u coloum vector
    end