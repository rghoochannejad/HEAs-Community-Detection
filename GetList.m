function list = GetList(nNode,A)

for i = 1:nNode
    R = find(A(i,:)==1);
    C = find(A(:,i)==1);
    
    list{i} = union(R,C);  
end


