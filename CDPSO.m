clc;
clear;
close all;

% ================== Community Detection Using  PSO  ===========
%% Load Graph
temp=load('totalAdj.mat');

%% initialization
A=temp.totalAdj;
nNode = size(A,1);
list = GetList(nNode,A);

c=1;  NeighborsNodes=[];
for i=1:size(A,1)
    for j=1:size(A,2)
         if A(i,j)==1
             m=1;
             NeighborsNodes(c,m)=i;
             m=m+1;
             NeighborsNodes(c,m)=j;
             c=c+1;
         end
    end
    
end

Problem.Edges = NeighborsNodes;
Problem.A = A;
Problem.nNode = nNode;
Problem.Position = 1:nNode;
%% Problem Definition
Problem.CostFunction = @(x) EvaluteGraph(x,Problem);      % Cost Function
Problem.List = list;


%% Set Params
Params.nPop = 200;
Params.MaxIt =150;

%% Community Detection Using EA
Results = CommunityDetectionUsingPSO(Problem,Params);


%% Results Evaluate
V = Results.F.Cluster;


%% Show Results
S = Problem.Position;
T = Results.F.GenoType;
g = digraph(S,T);

% G =Graph(g);
% figure(2)
% subplot(1,2,1),
plot(g,'NodeCData',V,'MarkerSize' ,5)
title([num2str(max(V)),' Clustering with PSO'])

save 'clusteringUsingPSO.mat'  'V';

