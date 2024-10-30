function [consensus] = compute_consensus(E, alpha, X)

f = E(X);
weight = exp(-alpha*(f-min(f)));
consensus = (1/sum(weight)).*sum((X.*weight),2);

end