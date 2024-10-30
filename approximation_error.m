clear; clc; close all;

global gamma g1 

lim1=4;

f1=@(x) 0;
g1=@(x) 1;

[X]=ndgrid(-lim1:0.01:lim1);
l=@(x)   (x.^2-2.2).^2-0.08.*x+0.5;
plot(X,l(X),'LineWidth', 2,'DisplayName','$f(x)$');
hold on;


Ntotal_set=[2,4,6,8];

for i=1:length(Ntotal_set)

 Ntotal=Ntotal_set(i);
 gamma=0.05;
 theta=0.5;

c=find_approximation(gamma,l,lim1,Ntotal,theta);
label_M = ['$V_n(x),M=$', num2str(Ntotal_set(i))];
approx=0;
for j=1:length(c)
    approx=approx+c(j)*phi1(X,j);
end

plot(X,approx,'LineWidth', 2, 'LineStyle', '-','DisplayName',label_M);
hold on;
end
legend('Interpreter','latex','FontSize',15)
xlabel('$x$','Interpreter','latex','FontSize',18) 
xlim([-3,3])
ylim([-1,6])

