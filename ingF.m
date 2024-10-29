function [wa]=ingF(i,d,Nf,NN,LELE,F0)
%global d Mf 
 wa=0;
            suma=0;
            for l=1:Nf
                pro=1;
                for m=1:d
                    %F0 is a tensor valued function from R^d--R^{d*n_f*d}
                    pro=pro*(LE(:,NN(i,m)+1)'*F0(:,m,l));
                end
                suma=suma+pro;
            end
            wa=wa+suma;