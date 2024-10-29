function [LE,LEP,LELE,LEPLE,LELELE,LEPLEPLE,LEPLELE]=legeval(gradmax)
global gp gw 

for i=1:gradmax+1
    LE(:,i)=leg(gp,i-1);
    LEP(:,i)=legp(gp,i-1);
end

LELE=zeros(length(gp),gradmax+1,gradmax+1);
LEPLE=zeros(length(gp),gradmax+1,gradmax+1);

for i=1:gradmax+1
 for j=1:gradmax+1
LELE(:,i,j)=gw'.*(LE(:,i).*LE(:,j));
LEPLE(:,i,j)=gw'.*(LEP(:,i).*LE(:,j));
 end
end


LELELE=zeros(gradmax+1,gradmax+1,gradmax+1);
LEPLELE=zeros(gradmax+1,gradmax+1,gradmax+1);
LEPLEPLE=zeros(gradmax+1,gradmax+1,gradmax+1);

for i=1:gradmax+1
 for j=1:gradmax+1
  for k=1:gradmax+1
LELELE(i,j,k)=gw*(LE(:,i).*LE(:,j).*LE(:,k));
LEPLELE(i,j,k)=gw*(LEP(:,i).*LE(:,j).*LE(:,k));
LEPLEPLE(i,j,k)=gw*(LEP(:,i).*LEP(:,j).*LE(:,k));
  end
 end
end