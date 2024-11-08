clear; clc; close all;

global Ntotal NN d W gamma Nf gp gw LE LEP LELE LEPLE F0 Ntotali LELELE LEPLELE LEPLEPLE lim

%% HJB solver

%Dimenion
d=2;

%Objective function:
objective_function='rastrigin';

if strcmp(objective_function,'ackley')
    E=@(v) -20*exp(-0.2*sqrt(1/d*sum(v.*v,1)))-exp(1/d*sum(cos(2*pi*v),1))+exp(1)+20;
    l=@(x,y) -20*exp(-0.2*sqrt(1/2.*(x.^2 + y.^2)))-exp(1/2.*(cos(2*pi*x)+cos(2*pi*y)))+exp(1)+20;
elseif strcmp(objective_function,'rastrigin')  
    E=@(v) sum(v.*v,1) + 10*sum(1-cos(2*pi*v),1);
    l=@(x,y) 20+x.^2 + y.^2 - 10*(cos(2*pi*x)+cos(2*pi*y));
end

%True minimizer
vstar=zeros(d,1); 

%epsilon: paramater for Tikhonov regularization
gamma=0.05;
epsilon=2*gamma;

%Gauss-Legendre quadrature
gp=[-0.960289856497536231684;-0.7966664774136267395916;-0.5255324099163289858177;-0.1834346424956498049395;0.1834346424956498049395;0.525532409916328985818;0.796666477413626739592;0.9602898564975362316836];
gw=[0.1012285362903762591525 0.222381034453374470544 0.313706645877887287338 0.3626837833783619829652 0.3626837833783619829652 0.31370664587788728734 0.222381034453374470544 0.1012285362903762591525];

%Region Omega [-lim,lim]^d
lim=2;

gp=lim*gp;
gw=lim*gw;

%The maimum degree in basis
gradmax=4;

%Generate the basis: full polynomial basis truncated by total degree M=gradmax
N_basis=gradmax*ones(1,d);
[Ntotal,NN,Ntotali]=baseinit(N_basis,gradmax,d);
[LE,LEP,LELE,LEPLE,LELELE,LEPLEPLE,LEPLELE]=legeval(gradmax);

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

%%  Controlled CBO

%Time horizon
T=10;
dt=0.01;
 
%Number of particles
N=[50,100,1000];

%Parameter for drift
lambda=1;
beta=1;

%Parameter for diffusion
sigma=0.7;
 
%Weight paramater for consensus
alpha=40;

%Paramater for initial distribution of particle system
upperbound=-0.5;
lowerbound=-1;

% Distance function
distance_control=NaN(length(N),1+T/dt);

% Variance function
Variance_control=NaN(length(N),1+T/dt);

for m = 1:length(N)

   % Initialization
    rng(1)
    X0=unifrnd(lowerbound,upperbound,d,N(m));
    X=X0;

    distance_control(m,1)=sum(vecnorm(X0-vstar).^2)/N(m);
    Expectation_control=sum(X0,2)/N(m);
    Variance_control(m,1)=1/2*sum(vecnorm(X0-Expectation_control).^2)/N(m);
   
    %Controlled CBO Algorithm
    for k = 1:T/dt
        %Compute consensus
        consensus=compute_consensus(E, alpha, X);
        %Controlled CBO update
        X=newcontrolledCBO_update(lambda,beta,dt,sigma,consensus,X,control,E,c);

        %distance function
        distance_control(m,k+1)=sum(vecnorm(X-vstar).^2)/N(m);
        % Variance
        Expectation_control=sum(X,2)/N(m);
        Variance_control(m,k+1)=1/2*sum(vecnorm(X-Expectation_control).^2)/N(m);
    end

end


%% The evolution of Var (\rho_t^N) and W_2(\rho_t^n ,\sigma_{x^*})for different number of particles N 

color_s=[0  0.4470  0.7410;0.8500  0.3250  0.0980;0.9290  0.6940  0.1250;0.4940  0.1840  0.5560;0.4660  0.6740  0.1880;0.3010  0.7450  0.9330;0.6350  0.0780  0.1840];
figure(1)
for m=1:length(N)      
    label_Var=['$\mathrm{Var}(\rho^N_t),N=$', num2str(N(m))];
    plot(0:dt:T,Variance_control(m,:),'LineWidth',4-m,'LineStyle','--',"color",color_s(m,:),'DisplayName',label_Var);
    hold on;
end

for m=1:length(N)
    label_distance=['$W_2^2(\rho^N_t,\delta_{x^*}),N=$', num2str(N(m))];
    plot(0:dt:T,distance_control(m,:),'LineWidth',4-m,'LineStyle','-',"color",color_s(m,:),'DisplayName',label_distance);
    hold on;
end

xlim([0,10])
set(gca,'YScale','log');

if strcmp(objective_function,'ackley')
ylim([10e-10,10e0])
elseif strcmp(objective_function,'rastrigin')  
ylim([10e-40,10e0])
end

ax.FontSize=13;
xlabel('$t$','Interpreter','latex','FontSize',15)
legend('Interpreter','latex','FontSize',18,'Location','northeast')
title('Controlled CBO')

%% Plot the trajectory of N = 50 particles
N=50;
figure(2)    

[X,Y]=meshgrid( -2:0.05:2,-2:0.05:2);
Z=reshape(l(X,Y), size(X));

if strcmp(objective_function,'ackley')
    s=surfc(X,Y,Z+1,'FaceAlpha',0.5);
    zlim([0,8]);
elseif strcmp(objective_function,'rastrigin')  
    s=surfc(X,Y,Z+10,'FaceAlpha',0.5);
    zlim([0,60]);
end

shading interp    % interpolate the colormap across the surface face
colorbar
hold on;

[X_grid,Y_grid]=ndgrid(-2:0.05:2,-2:0.05:2);
contourf(X_grid,Y_grid,l(X_grid,Y_grid),'FaceAlpha',0.5,'LineStyle','none');
colorbar;
hold on;

view(-20,24)
lightangle(-20,25)
set(gca,'FontSize',13);

X_plot=NaN(N,2,T/dt+1);
rng(1)
X0=unifrnd(lowerbound,upperbound,d,N);

 for i=1:N
    plot(X0(1,i), X0(2,i),'o-','MarkerFaceColor','b','MarkerSize',7)
    hold on;
    X=X0;
    %Controlled CBO
    for k = 1:T/dt
        X_plot(i,:,k)=X(:,i);
        %Compute consensus
        consensus=compute_consensus(E, alpha, X);
        %Controlled CBO update 
        X=newcontrolledCBO_update(lambda,beta,dt,sigma,consensus,X,control,E,c);
    end
     w=rem(i,6)+1;
     plot(squeeze(X_plot(i,1,:)),squeeze(X_plot(i,2,:)),'o-','MarkerSize',2,'MarkerFaceColor',color_s(w,:))
     plot(squeeze(X_plot(i,1,T/dt)),squeeze(X_plot(i,2,T/dt)),'-o','MarkerFaceColor','red','MarkerSize',8)
end


if upperbound==0.5
    xlim([-1.5,1]);
    ylim([-1.5,1]);
elseif upperbound==-0.5
   xlim([-1.5,0.5]);
   ylim([-1.5,0.5]);
end
title('Trajectory of controlled CBO')

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