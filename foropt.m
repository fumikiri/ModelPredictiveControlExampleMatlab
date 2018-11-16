function [bound_ofUmax,bound_ofUmin,bound_ofYmax,bound_ofYmin,H,M1,M3,C1] = foropt(phi,phi_phi,R,Umax,Umin,Ymax,Ymin,in,out,Np,Nc)
if size(Umax)~=0
bound_ofUmin=[];
bound_ofUmax=[];
for j=1:Nc
   bound_ofUmin((j-1)*in+1:j*in,:)=Umin;%[-40.68;-2.102;0];%[-500;-500;-500]
   bound_ofUmax((j-1)*in+1:j*in,:)=Umax;%[79.32;4.898;10];%[500;500;500]
end
C1=zeros(Nc*in,in);
C2=zeros(Nc*in,Nc*in);
      for j=1:Nc;
          C1(in*(j-1)+1:in*j,:)=eye(in,in);
      end
      for j=1:Nc
          C2(in*(j-1)+1:in*Nc,in*(j-1)+1:in*j)=C1(1:in*(Nc-j+1),:);
      end
M1=[-C2;C2];
else
    M1=[];
end
if size(Ymax)~=0
    for j=1:Np
   bound_ofYmin((j-1)*out+1:j*out,:)=Ymin;%[-40.68;-2.102;0];%[-500;-500;-500]
   bound_ofYmax((j-1)*out+1:j*out,:)=Ymax;%[79.32;4.898;10];%[500;500;500]
    end
    M3=[-phi;phi];
else
    M3=[];
end
H=2*(phi_phi+R*eye(Nc*in,Nc*in));
end

