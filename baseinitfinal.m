function [Ntotal1,NN1]=baseinitfinal(gradmax,even,basisgen)
global d

switch basisgen

%% total degree
case 1 

B=sparse([zeros(1,d);speye(d)]);
B=B';
grad=1;
while(grad<=gradmax)
    A=B;
    %B=sparse(1,d);
    nn=size(A,2);
    tic
for i=1:nn
    well=sparse(d,nn);
    for j=1:nn
        wa=A(:,i)+A(:,j);
        if (sum(wa)==grad)
        well(:,j)=wa;
        end
    end
    wellt=unique(well','rows');
    B=[B wellt'];
end
toc
B=unique(B','rows');
B=B';
grad=grad+1

end

NN1=B';
%NN1(1,:)=[];
Ntotal1=size(NN1,1);


grados=sum(NN1,2);
lista=zeros(Ntotal1,1);

if (even==1)
for i=1:Ntotal1
    par=mod(grados(i),2);
    if (par==1 ||grados(i)>gradmax) 
    %if (grados(i)==3 || grados(i)>gradmax) 
    %if (grados(i)>gradmax || grados(i)<=1) 
        %lista=[lista i];
        lista(i)=1;
    end
end

else
    for i=1:Ntotal1
    par=mod(grados(i),2);
    %if (par==1 ||grados(i)>gradmax) 
    %if (grados(i)==3 || grados(i)>gradmax) 
    if (grados(i)>gradmax || grados(i)<1) 
        %lista=[lista i];
        lista(i)=1;
    end
    end
    
end
NN1(lista==1,:)=[];
Ntotal1=length(NN1);

%% Hyberpolic Cross Fast

case 2
N1=Smolyak_Elem_Isotrop(d,gradmax);
gradpol=2^gradmax;
%gradpol=gradmax;
grados=prod(N1,2);
N1(grados>gradpol+1,:)=[];
NN1=N1-1;
Ntotal1=size(NN1,1);


%% Euclidean degree
case 3
B=sparse([zeros(1,d);speye(d)]);
B=B';
grad=1;
while(grad<=gradmax)
    A=B;
    %B=sparse(1,d);
    nn=size(A,2);
    tic
for i=1:nn
    well=sparse(d,nn);
    for j=1:nn
        wa=A(:,i)+A(:,j);
        if (sum(wa)==grad)
        well(:,j)=wa;
        end
    end
    wellt=unique(well','rows');
    B=[B wellt'];
end
toc
B=unique(B','rows');
B=B';
grad=grad+1

end

NN1=B';
%NN1(1,:)=[];
Ntotal1=size(NN1,1);

%%%euclidean degree%%%
lista=zeros(Ntotal1,1);
gradosh=sqrt(sum(NN1.^2,2));

for i=1:Ntotal1
    if gradosh(i)>(gradmax)
        lista(i)=1;
    end
end
NN1(lista==1,:)=[];
Ntotal1=length(NN1);

%% Smolyak
case 4      
NN1=Smolyak_Elem_Isotrop(d,gradmax)-1;
Ntotal1=size(NN1,1);

end

