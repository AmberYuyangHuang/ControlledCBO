function [wa]=inG(i,j,Nu,LELE,LEPLE,gp,d,NN)
wa=0;

        for p=1:d
            for l=1:Nu
                prod=1;
                for m=1:d
                switch m
                    case p
                        fact=LEPLE(:,NN(i,p)+1,NN(j,p)+1)'*zeros(size(gp));
                    otherwise
                        fact=LELE(:,NN(i,m)+1,NN(j,m)+1)'*u0sep(l,m,gp);
                end
                prod=prod*fact;
                end
                wa=wa+prod;
            end
        end