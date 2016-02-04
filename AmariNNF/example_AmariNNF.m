%% Defaults stumulus, functions and their parameters
A=3;
B=2;
C=0.6;
a=1;
b=3;
d=4;
h=1;

t =  0:0.1:200;
x = -8:0.1:8; Mx = length(x); % default discretization used in the code

w = ( A*( abs(x)<=a ) - B*( ( (abs(x) > a) + (abs(x) <= b) ) ==2) );
wu = 0.8 * x .* (abs(x)<=b);
Stim = C*(1-abs(x)/d).* ( abs(x) <=d);

Uin = -h * ones(Mx,1);

%% Simulation
[Uxt1] = AmariNNF(w,Uin,Stim); % [Uxt]=AmariNNF(C,w,ini,Stim)
[Uxt2] = AmariNNF(w,randn(Mx,1),Stim);
h=figure
subplot(2,2,3),
     [tt,xx] =  meshgrid(t,x);
     mesh(tt,xx,Uxt1)
     xlabel('time')
     ylabel('location')    
     title('U(x,t), with U_0=-h') % for (1.6)
subplot(2,2,4 ),
     [tt,xx] =  meshgrid(t,x);
     mesh(tt,xx,Uxt2)
     xlabel('time')
     ylabel('location')    
     title('U(x,t), with U_0=GWN') % for (1.6)
         
subplot(2,2,[1 2]), 
        plot(x,Uxt1(:,end),'-or',x,Uxt2(:,end),'-g',x,zeros(1,Mx)), title('Ex2.1-1')
        legend('U_0=-h','U_0=GWN','zeros level'), xlabel('x'), ylabel('U(x,infinity)')
        title('U(x,infinity')