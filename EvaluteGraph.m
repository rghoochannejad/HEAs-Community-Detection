function [z,Out] = EvaluteGraph(x,Data)
A = Data.A;
Position = Data.Position;
% GG = Data.V0;

X = [Position',x'];
sg = graph(X(:,1),X(:,2));
% figure,plot(sg)
V = conncomp(sg);
nc = numel(unique(V));
if nc>1
    z =QFModul(V,A);
   
else
    z = 1e5;
end
    Out.z = z;
    
    Out.GenoType = x;
    Out.Cluster = V;
end

