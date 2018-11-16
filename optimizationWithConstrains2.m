%function [] = optimizationWithConstrains(phi,phi_phi,R,Umax,Umin,Ymax,Ymin,in,out,Np,Nc)
function [M2, N2] = optimizationWithConstrains(deltaumax,deltaumin,Np,Nc)

M2 = [-eye(Nc)
       eye(Nc)];
N2 = [-ones(Nc,1)*deltaumin
      ones(Nc,1)*deltaumax];

























