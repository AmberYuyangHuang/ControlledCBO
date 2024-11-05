clear; clc; close all;

global Ntotal NN d W gamma Nf gp gw LE LEP LELE LEPLE F0 Ntotali LELELE LEPLELE LEPLEPLE lim

%% HJB solver

%Dimenion
d=2;

%Objective function:
objective_function='rastrigin';

if strcmp(objective_function,'ackley')
    E=@(v) -20*exp(-0.2*sqrt((1/d).*sum(v.*v,1)))-exp((1/d).*sum(cos(2*pi*v),1))+exp(1)+20;
elseif strcmp(objective_function,'rastrigin')    
    E=@(v) sum(v.*v,1) + 10*sum(1-cos(2*pi*v),1);
end

%True minimizer
vstar=zeros(d,1); 

%parameters of HJB solver
gamma=0.05;
epsilon=2*gamma;

%Gauss-Legendre quadrature
gp=[-0.960289856497536231684;-0.7966664774136267395916;-0.5255324099163289858177;-0.1834346424956498049395;0.1834346424956498049395;0.525532409916328985818;0.796666477413626739592;0.9602898564975362316836];
gw=[0.1012285362903762591525 0.222381034453374470544 0.313706645877887287338 0.3626837833783619829652 0.3626837833783619829652 0.31370664587788728734 0.222381034453374470544 0.1012285362903762591525];

%Region Omega [-lim,lim]^d
lim=2;

gp=lim*gp;
gw=lim*gw;

%Maximum degree in basis 
gradmax= 2;

N_basis=gradmax*ones(1,d);

%Generate the basis
basis_type='hyperbolic_cross';
if strcmp(basis_type,'hyperbolic_cross')
    [Ntotal,NN]=baseinitfinal(gradmax,0,2);
    NN(1,:) = [];
    Ntotal=Ntotal-1;
    Ntotali=polycomb(1:Ntotal,1:Ntotal)';
    [LE,LEP,LELE,LEPLE,LELELE,LEPLEPLE,LEPLELE]=legeval(2^gradmax);
elseif strcmp(basis_type,'full')
    %[Ntotal,NN]=baseinitfinal(gradmax,0,1);
    [Ntotal,NN,Ntotali]=baseinit(N_basis,gradmax,d);
    [LE,LEP,LELE,LEPLE,LELELE,LEPLEPLE,LEPLELE]=legeval(gradmax);
end

Megawind=polycomb(1:Ntotal,1:Ntotal,1:Ntotal)';
Megadim=Ntotal*Ntotal*Ntotal;
[W]=Weg(Megawind,Megadim,LELELE,LEPLEPLE);

%Coefficients of polynomial approximation V_n
if strcmp(objective_function,'rastrigin')
    [F0,Nf]=f0parab();
    %if the objective function is separable
    solution=approximation(E);
    a=cost_new_t();
elseif strcmp(objective_function,'ackley')  
    %if the objective function is non-separable
    solution=approximation_ns(E);
    a=cost_new_mc(E);
end

%Control function
control=@(x) compute_control(x,epsilon,solution);

%Coefficients of f^approx
c=zeros(Ntotal,1);
b=bb();

for i=1:length(Ntotal)
    c(i)=a(i)/b(i);
end

%% Controlled CBO

%Time horizon
T=10;
dt=0.01;
 
%Number of particles
N=50;

upperbound=-0.5;
lowerbound=-1;

%Parameter for drift
lambda=1;
beta=1;

%Parameter for diffusion
sigma=0.7;
 
%Weight paramater for consensus
alpha=40;

total_num=100;
sum_distance=0;

for num=1:total_num
    % Initialization of particle system
    X0=unifrnd(lowerbound,upperbound,d,N);
    X=X0;
   
    % Controlled CBO Algorithm
    for k = 1:T/dt
        %Compute consensus 
        consensus=compute_consensus(E, alpha, X);
        %Controlled CBO update
        X=newcontrolledCBO_update(lambda,beta,dt,sigma,consensus,X,control,E,c);    
    end
    sum_distance=sum_distance+sum(vecnorm(X-vstar).^2)/N;
end

fprintf('Averaged distance: %e\n', sum_distance/total_num)



%% Controlled CBO update

function [X] = newcontrolledCBO_update(lambda,beta,dt,sigma,consensus,X,control,E,c)

global NN d Ntotal

N=size(X,2);

approx=zeros(1,N);

for i=1:N
    phi=zeros(d,Ntotal);
    for j=1:Ntotal
        prod=1;
        for p=1:d
            prod=prod*leg(X(p,i),NN(j,p));
        end
        phi(i,j)=prod;
    end
    %Compute f_approximated
    approx(1,i)=0;
    for q=1:Ntotal
        approx(1,i)=approx(1,i)+c(q)*phi(i,q);
    end
end

%Drift term
    for i =1:N
        if  E(consensus)<E(X(:,i))
            X(:,i)=X(:,i)-lambda*(X(:,i)-consensus)*dt;
        end
       
        if  approx(i)<E(X(:,i))
            X(:,i)=X(:,i)+beta*control(X(:,i))*dt;
        end
    end

%Diffusion term
X=X+sigma*abs(X-consensus*ones(1,N))*sqrt(dt).*randn(d,N);

end
