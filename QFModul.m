function Q = QFModul(V,A)

% Modularity  quality function

m = sum(sum(A));
Q = 0;
COMu = unique(V);
for j=1:length(COMu)
    Cj = find(V==COMu(j));
    lc = sum(sum(A(Cj,Cj)));
    dc = sum(sum(A(Cj,:)));
    if dc>0
        Q = Q + lc/m-(dc/(2*m))^2;
    end
end
Q = -Q;

