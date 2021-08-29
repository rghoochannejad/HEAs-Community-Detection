function [gBest,MPS] = GetGBest(MPS)

[n,m] = size(MPS);
for i = 1:n
    
    Li = MPS(i,:);
    j = 0;
    while (1)
        j = j + 1;
        if j<size(MPS,1)
            if ~(j==i)
                Lj = MPS(j,:);
                
                if sum(Li==Lj)==m
                    MPS(j,:) = [];
                    j = j -1;
                end
            end
        else
            break;
        end
        
    end 
    if i>=size(MPS,1)
        break;
    end
end

MP0 = MPS(1,:)';
n = size(MPS,1);
if n==1
    gBest = MP0';
    return;
end
for i = 2:n
    MPi = MPS(i,:)';
    wi = pinv(MPi'*MPi)*(MPi'*MP0);
    Vi = MPi*wi;
    MP0 = ((i - 1)/i)*MP0 + (1/i)*Vi;
end

gBest = MP0';

end



