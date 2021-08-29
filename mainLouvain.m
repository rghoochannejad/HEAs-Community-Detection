clc;
clear;
close all;

%% Calculate Jaccard similarity

X = xlsread('CosineSimilarity.xls','sheet2');

X=X(2:end,2:end);

for i= 1:size(X,1)
     a = X(i,:);
   for j= 1:size(X,1)
  if i==j
      continue;
  else
      b = X(j,:);
    JD(i,j) =sum(a & b)/sum(a | b) ;  %jaccard similarity
    JD(isnan(JD))=0 ;  
%    w JI = 1 - JD;  % jaccard distance 
  end
   end
end
% a = X(1,:);
% b = X(2,:);
% JD =sum(a & b)/sum(a | b) ;
%  JD = pdist(X,'jaccard') ;  % jaccard distance
%  JI = 1 - JD;   

%% Dataset Normilization 
read = xlsread('data.xlsx') ;
Data=read(:,:); 
minData = min(Data(:));
maxData = max(Data(:));
scaled  = (Data - minData) / (maxData - minData);  % Scaled to [0, 1]

%% Calculate Cosine similarity

n_row = size(Data,1);
norm_r = sqrt(sum(abs(Data).^2,2)); 
S2 = zeros(n_row,n_row);
%  dot product is a matrix product of a row vector with a column vector
for i = 1:n_row
  for j = i:n_row
    S2(i,j) = dot(Data(i,:), Data(j,:)) / (norm_r(i) * norm_r(j));
    S2(j,i) = S2(i,j);
  end
end
%  xlswrite('CosineSimilarity_new.xls',S2) ;

%%  Calculate CloseNess Centrality

count=1;
n=size(X,1) ;
for i=1:n
    for j=1:n
        if(X(i,j) ~=0)
       m1(count,1)=i;
       m2(count,1)=j;
       weights(count,1)=X(i,j);
       count=count+1;
        end
    end
end


% G=graph(m1,m2,weights);
% S=G.Edges;
% save 'Tempadjacency.mat' 'S'

% closeNess=centrality(G,'closeness');
% BNess=centrality(G,'betweenness');
% 
% xlswrite('closeNess.xls',closeNess);
% xlswrite('betweenNess.xls',BNess);

% plot(G,'XData',x,'YData',y,'ZData',z,'EdgeLabel',G.Edges.Weight)
% LWidths = 5*G.Edges.Weight/max(G.Edges.Weight);
% plot(G,'EdgeLabel',G.Edges.Weight,'LineWidth',LWidths)
% saveas(gcf,'graph.jpg')

%% Community detection with Louvain method
myData=xlsread('CosineSimilarity.xls','sheet1') ;
myData=myData(2:end,2:end);
alpha=0.7;
myMatrix=alpha .* JD +(1-alpha).* myData ;
%% =================== Set a threshold value ====================
Threshold= input('Enter a Threshold');
for i=1:size(myMatrix,1)
    for j=1:size(myMatrix,2)
        if myMatrix >=Threshold
        myMatrix(i,j)=myMatrix(i,j);
        else
           myMatrix(i,j)=0; 
        end
    end
end
%% =================== Creating a graph ==============
count=1;
n1=size(myMatrix,1);
for i=1:n1
    for j=1:n1
        if(myMatrix(i,j) ~=0)
       v1(count,1)=i;
       v2(count,1)=j;
       weights(count,1)=myMatrix(i,j);
       count=count+1;
        end
    end
end

 G2=graph(v1,v2,weights);

% adj = adjacency(G);
% totalAdj = full(adj);


% S=load ('Tempadjacency');
% S = table2array(S.S);
% G = graph(S(:,1), S(:,2));
% plot(G);
% adj = adjacency(G);

% gamma = 1;
% k = full(sum(totalAdj));
% twom = sum(k);
% B = full(totalAdj - gamma*k'*k/twom);
% [S,Q] = genlouvain(B);   % S means each Node is assigned to one cluster
% save 'S.mat' 'S'
% Q = Q/twom               % quality of modularity
 %============ showing each node into itself cluster =====================
%%
 clc;
clear;
close all;

S=load('S.mat');
S=S.S;
numCluster=max(S);
count=1;
for i=1:numCluster
    
  for j=1:size(S,1)
     if(i==S(j))
         clusters(i,count)=j;
         count=count+1;
     end
     
  end
  count=1;
end
% plot communities
PlotCommunities(clusters,numCluster);
