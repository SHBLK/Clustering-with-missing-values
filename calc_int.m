function [logprob] = calc_int(G,g,mg,Sg,ng,v,m,Sigma)

d = length(m);
logdet = 0;
trc = 0;
A = v*inv(Sigma);
b = v*(Sigma\m);
q = m'*b;
for gg=1:G
    idx = g(gg,:)==1;
    logdet = logdet + ng(gg)*log(max(det(2*pi*Sigma(idx,idx)),realmin));
    trc = trc + trace(Sg{gg}/Sigma(idx,idx));
    temp = inv(Sigma(idx,idx));
    q = q + ng(gg)*mg{gg}'*temp*mg{gg};
    tt = zeros(d);
    tt(idx,idx) = temp;
    A = A + ng(gg)*tt;
    tt = zeros(d,sum(idx));
    tt(idx,:) = temp;
    b = b + ng(gg)*tt*mg{gg};
end
q = q - b'*inv(A)*b;
logprob = -0.5*(-d*log(v)+log(max(det(2*pi*Sigma),realmin))...
        + log(max(det(A/(2*pi)),realmin))+logdet+trc+q);

end

