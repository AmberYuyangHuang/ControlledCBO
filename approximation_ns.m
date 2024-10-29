%close all
%clear all
%clc
function [sol]=approximation_ns(E)

global Ntotal NN d W gamma LELE


%discount factor
lambda_ad=0.1;
epsilon=2*gamma;

%Nu=2*d;
%[U0]=uinit(Nu);
%[G0]=grad0(Ntotali,Ntotal,NN,Nu,LELE,LEPLE,gp,d);

[U0]=zeros(Ntotal,1);
[G0]=zeros(Ntotal,Ntotal);

[L]=cost_new_mc(E);
[M]=mass(Ntotal,NN,d,LELE);

A=-lambda_ad*M+G0;
b=-(L+U0);
A=A.*(abs(A)>1e-10);
sol=A\b;
c=0.*sol;
%lambda_ad=0.00002;
while(lambda_ad>0.00001)
while(norm(sol-c)>1e-8)
c=sol;
[G]=gradit(c);
[U1]=uti(c);
A=-lambda_ad*M+G;
b=-(L+U1);
A=A.*(abs(A)>1e-10);
sol=A\b;
norm(sol-c);
end
lambda_ad=0.5*lambda_ad;
end
%sol;
end
