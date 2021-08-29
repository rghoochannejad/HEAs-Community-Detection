function Results = CommunityDetectionUsingPSO(Problem,Params)

%% Problem Definition

CostFunction=Problem.CostFunction;        % Cost Function

nVar=Problem.nNode;             % Number of Decision Variables

VarSize=[1 nVar];   % Size of Decision Variables Matrix

VarMin= 1;          % Lower Bound of Variables
VarMax= Problem.nNode;          % Upper Bound of Variables

List = Problem.List;
%% PSO Parameters

MaxIt=Params.MaxIt;      % Maximum Number of Iterations

nPop=Params.nPop;        % Population Size (Swarm Size)


% c1 = 1.496;        % Personal Learning Coefficient
c1=1.5 ;
c2 =1.5;        % Global Learning Coefficient
roh = 0.75;
Tmax = 10;
MPS = [];
wmin = 0.1;
wmax = 1;
% Velocity Limits
VelMax=0.1*(VarMax-VarMin);
VelMin=-VelMax;
tic

%% Initialization

empty_particle.Position=[];
empty_particle.Cost=inf;
empty_particle.Velocity=[];
empty_particle.Best.Position=[];
empty_particle.Best.Cost=[];

particle=repmat(empty_particle,nPop,1);

BestSol.Cost=inf;

GlobalBest = empty_particle;

for i=1:nPop
    % Initialize Position
    NewPosition = InitialPop(Problem);
    particle(i).Velocity = unifrnd(.2*VarMin,.2*VarMax,VarSize);
    % Evaluation
    Cost=CostFunction(NewPosition);
    
    % Update Personal Best
    if Cost<particle(i).Cost
        particle(i).Cost = Cost;
        particle(i).Position = NewPosition;
        particle(i).Best.Cost  = Cost;
        particle(i).Best.Position = NewPosition;
    end
    % Update Global Best
    if  Cost<GlobalBest.Cost
        GlobalBest.Position = NewPosition;
        GlobalBest.Cost= Cost;
        BC = Cost;
        gBest = GlobalBest.Position;
    end
end


NFE = nPop; % Number of Function Evaluation
BestCost = zeros(MaxIt,1);
nfe = zeros(MaxIt,1);
MaxRun=10;
Results=[];
%% ===================== PSO Main Algorithm  ======================
for run=1:MaxRun
 disp (['============ run( ' num2str(run) ' ) ================ ' ])
for it=1:MaxIt
    w=(wmax - wmin)*(1 - it/MaxIt) + wmin;
    for i=1:nPop
        
        % Update Velocity
        particle(i).Velocity = w*particle(i).Velocity ...
            +c1*rand(VarSize).*(particle(i).Best.Position-particle(i).Position) ...
            +c2*rand(VarSize).*(gBest-particle(i).Position);
        
        expV = exp(-particle(i).Velocity);
        signV = abs((1- expV)./(1 + expV));
        Position = particle(i).Position;
        BP = GlobalBest.Position;
        for v = 1:nVar
            K = List{v};
            nj = numel(K);
            if .5<signV(v) && nj>1
                
                k = randsample(nj,1);
                Position(v) = K(k) ;
                
            else
%                Position(v) = BP(v);
            end
        end
%         [Position,flag] = CheckVisibleEdge(List,Position);
        % Update Position
        particle(i).Position = Position;
        
        % Velocity Mirror Effect
        % Apply Velocity Limits
        particle(i).Velocity = min(max(particle(i).Velocity,VelMin),VelMax);
        
        
        % Evaluation
        particle(i).Cost = CostFunction(particle(i).Position);
        
        % Update Personal Best
        if particle(i).Cost<particle(i).Best.Cost
            
            particle(i).Best.Position=particle(i).Position;
            particle(i).Best.Cost=particle(i).Cost;
            
            % Update Global Best
            if particle(i).Best.Cost<BestSol.Cost
                
                BestSol=particle(i).Best;
                
            end
            
            
        end
        
    end
    
    % Update Global Best
    if BC>BestSol.Cost
        gBest = BestSol.Position;
        MPS = [];Tmax = 0;
    elseif BC == BestSol.Cost
        Tmax = Tmax + 1;
        MPS = [MPS;BestSol.Position];
    elseif Tmax>2
        MPS = [MPS;particle(randi(nVar,1)).Position];
    end
    
    if (Tmax>2)
        [gBest,MPS] = GetGBest(MPS);
    end
    
    BestCost(it)=BestSol.Cost;
    BC = BestCost(it);
    disp(['Iteration ' num2str(it) ': Best Cost = ' num2str(BestCost(it))]);
        
end
array(run,:) = GlobalBest.Position;

end %------ end of run

disp('End of PSO...!');
ResultMean = mean(array);
Time= toc;
 Results.Time = Time;
 Results.NFE = nfe(end);
p = ResultMean;
[~,Out] = CostFunction(p);
Results.F = Out;
 % masalan dar row aval onai ke meghdar 1 darand yani oon node ha
 % motealegh be cluster 1 hastand.

end