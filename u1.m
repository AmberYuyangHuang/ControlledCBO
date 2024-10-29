function [val]=u1(x,c,gamma)
global g1
val=0;
for i=1:length(c)
    val=val+c(i)*phipx1(x,i);
end
val=-(0.5/gamma)*g1(x)*val;
    
