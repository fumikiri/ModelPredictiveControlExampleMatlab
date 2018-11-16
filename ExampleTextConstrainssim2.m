clear 
clc
close all

% E = [2 -1
%      -1 1];
% F = [-1
%      0]
%  
% M = [-1  0
%       0 -1
%       3  2]
% gamma = [0 
%          0
%          4]
% x0 = -E\F
% 
% 
% eta=QPhild(E,F,M,gamma)





Ts=2e-1; % Sampling time 1ms
Tsim=50; % Simulation time 20s
t=0:Ts:Tsim-1*Ts;
Nsim=Tsim/Ts;



% Plant
plant = tf(10,[1 0.1 3]);
[A B C D]=tf2ss(plant.num{1},plant.den{1})
sys=ss(A,B,C,D);
plantd = c2d(sys,Ts)
xplant(1:2,1:Nsim)=0;
yplant(1,1:Nsim)=0;

[out,state] = size(plantd.C);
[~,in] = size(plantd.B);


Np = 25;
Nc = 10;
Rs = 1;
Rw = 4e-2*1;
Q= 1
Rsbar=ones(Np,1);
[F,phi,phi_phi,phi_f,phi_r] = mpcgainOrigin(plantd.A,plantd.B,plantd.C,Nc,Np,Rs,Q)
[nnn mmm]=size(phi_phi)
Umax = 1;
Umin = -1;
Ymax = 1.3;
Ymin = -0.1;
deltaumax = 0.5;
deltaumin = -1;

[bound_ofUmax,bound_ofUmin,bound_ofYmax,bound_ofYmin,H,M1,M3,C1] = foropt(phi,phi_phi,Rw,Umax,Umin,Ymax,Ymin,in,out,Np,Nc);


ref(1:Nsim/2)   = 0;
ref(Nsim/4+1:Nsim/4*3) = 1;
ref(Nsim/4*3+1:Nsim) = 0;

u(1:Nsim) = 0;
deltau(1:Nsim) = 0;


% Simulation
for i=2:Nsim-1
    Xh = [(xplant(:,i)-xplant(:,i-1));yplant(1,i-1)];
    %deltaU=((phi_phi+eye(nnn,mmm)*Rw)\phi'*(Rsbar*ref(i)-F*Xh));
    
    
    [M2, N2] = optimizationWithConstrains2(deltaumax,deltaumin,Np,Nc);
    A_cons = [M2];
    b = [N2];
    H = 2*(phi_phi+eye(nnn,mmm)*Rw);
    f = -2*phi'*(Rsbar*ref(i)-F*Xh);
    eta=QPhild(H,f,A_cons,b);
    deltau(i)=([1 zeros(1,Nc-1)])*eta;
    u(i) = u(i-1)+deltau(i);
%     
     %deltau=([1 zeros(1,Nc-1)])*((phi_phi+eye(nnn,mmm)*Rw)\phi'*(Rsbar*ref(i)-F*Xh));
    
    %u(i)=u(i-1)+deltau;
    
    %u(i) = ref(i);
    xplant(:,i+1)=plantd.A*xplant(:,i)+plantd.B*u(i);
    yplant(1,i)=(plantd.C*xplant(:,i)+plantd.D*u(i));
end


figure;
subplot(3,1,1)
hold on
plot(t,yplant(1,:))
plot(t,ref(1,:))
subplot(3,1,2)
plot(t,u)
subplot(3,1,3)
plot(t,deltau)




