
function [c]=find_approximation(gamma,l,lim1,Ntotal,theta)


lambda_ad=theta;

f1=@(x) 0;
g1=@(x) 1;

M=sparse(Ntotal,Ntotal);

for i=1:Ntotal
    for j=1:Ntotal
    M(i,j)=integral(@(x) phi1(x,i).*phi1(x,j),-lim1,lim1);
    end
end

K1=sparse(Ntotal,Ntotal);
K20=sparse(Ntotal,Ntotal);
F=zeros(Ntotal,1);
u01=0;
for i=1:Ntotal
    for j=1:Ntotal
        K20(i,j)=integral(@(x) phi1(x,i).*u01.*(phipx1(x,j).*g1(x)),-lim1,lim1);

    end
end

for i=1:Ntotal
F(i)=-integral(@(x) phi1(x,i).*(l(x)+gamma*u01.^2),-lim1,lim1);
end

A=-lambda_ad*M+(K1+K20);
A=A.*(abs(A)>1e-8);
sol=A\F;

c=sol;

K2=zeros(Ntotal,Ntotal,Ntotal);
K2i=zeros(Ntotal,Ntotal);
for i=1:Ntotal
    for j=1:Ntotal
        for k=1:Ntotal
            K2(i,j,k)=(-0.5/gamma)*integral(@(x) phi1(x,i).*(phipx1(x,j).*g1(x)).*(phipx1(x,k).*g1(x)),-lim1,lim1);
        end
    end
end

it=0;
while(lambda_ad>0.00001)
error=1;
while(error>1e-10)

K2i=0*K2i;
for i=1:Ntotal
        K2i=K2i+squeeze(K2(:,:,i))*c(i);
end


for i=1:Ntotal
F(i)=-integral(@(x) phi1(x,i).*(l(x)+gamma*u1(x,c,gamma).^2),-lim1,lim1);
end
A=-lambda_ad*M+(K1+K2i);

sol=A\F;

error=norm(sol-c);
c=sol;
it=it+1;
end
lambda_ad=lambda_ad.*0.5;
end


end