function [Uxt]=AmariNNF(w,Uin,Stim)
%% Exercise 8_2: Simulate a Nonlinear Neural Field
% C for amp of stimulus, 
% w for interaction kernel,
% ini for intitial condition over space
% Stim for stimulus over space

% close all
% clear 
%% Parameters
tau = 10;
% discrete steps
T=200;
dx = 0.1; 
x = -8:dx:8;
Mx = length(x);
dt = 0.1;
t = 0:dt:T;
Nt = length(t);

Uxt = zeros(Mx,Nt);% Solution
Uxt(:,1) = Uin; % Put in initial Condition

%%  Simulation of the NF: u' = -u + conv_x (w * I(u)) + s(x) - h
for t_step=2:Nt
    Uxt_told = Uxt(:,t_step-1);
    Uxt_step =  ( sign(Uxt_told) .* ( Uxt_told>0 ) ); % Riemann Sum approximation: int w(x-x') * I(u(x',t)) dx'                
    
    % **Option1** for convolution, directly computing convolution of all possible locations
    convWU = dx * conv(w,Uxt_step,'same'); % cf: https://www.youtube.com/watch?v=LZ0qjZezGkQ    
    % **Option2** for convolution
    for x_step=1:Mx                
        X = x(x_step) - x; % X = x-x', in the integral term        
        % W = wsca * ( A*( abs(X)<=a ) - B*( ( (abs(X) > a) + (abs(X) <= b) ) ==2) ); % interaction Kernel w
        % W2 = ( 0.8 * X.^2 .* (abs(X)<=b) );
        % plot(W) % Kernel (1.6) shape exploration                                
        %% FWD Euler Approxi
        Uxt(x_step,t_step) =  Uxt_told(x_step)...                                                      
                            + (dt/tau) * ( -Uxt_told(x_step) + Stim(x_step) + convWU(x_step) -1 );                                                     
                              
    end
end    