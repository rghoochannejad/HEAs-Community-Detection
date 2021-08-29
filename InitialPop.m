function p = InitialPop(Data)
n = Data.nNode;
List = Data.List;

p = zeros(1,n);
for i = 1:n
    K = List{i};
    if(isempty(K))
       
      p(i)=i;  
    else
    nj = numel(K);
    k = randsample(nj,1);
    p(i) = K(k);
    end
end


%  [p,flag] = CheckVisibleEdge(List,p);
end