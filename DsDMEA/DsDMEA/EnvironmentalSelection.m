function Population = EnvironmentalSelection(Population,W,W1,N)

    %% Non-dominated sorting
    %[FrontNo,MaxFNo] = NDSort(Population.objs,N);
    %Next = FrontNo < MaxFNo;
%     Fitness = CalFitnessobjobj(Population.objs);
    Fitness = CalFitnessnewdec(Population.objs,Population.decs);
    Next = Fitness<1;
    if sum(Next)<N
        [~,Rank] = sort(Fitness);
        Next(Rank(1:N)) = true;
    else
    
    %% Select the solutions in the last front
    Last   = find(Next);
    Choose = LastSelection(Population(Last).decs,Population(Last).objs,W,W1,N);
    Next(Last(~Choose)) = false;
    end
    % Population for next generation
    Population = Population(Next);
end

function Remain = LastSelection(PopDec,PopObj,W,W1,K)
% Select part of the solutions in the last front

    N  = size(PopDec,1);
    %NW = size(W,1);
    NW = length(W);
    Wobj = W.objs;
    Wdec = W1.decs;
    %{
    %% Normalization
    %% 决策空间参考点归一化
    Dmin = min(W1.decs,[],1);
    Dmax = max(W1.decs,[],1);
    Wdec=(W1.decs-repmat(Dmin,NW,1))./repmat(Dmax-Dmin,NW,1);
    %% 目标空间参考点归一化
    Omin = min(W.objs,[],1);
    Omax = max(W.objs,[],1);
    Wobj =(W.objs-repmat(Omin,NW,1))./repmat(Omax-Omin,NW,1);
    %% 决策变量归一化
    PDmin = min(PopDec,[],1);
    PDmax = max(PopDec,[],1);
    PopDec=(PopDec-repmat(PDmin,N,1))./repmat(PDmax-PDmin,N,1);
    %% 目标向量归一化
    POmin = min(PopObj,[],1);
    POmax = max(PopObj,[],1);
    PopObj =(PopObj-repmat(POmin,N,1))./repmat(POmax-POmin,N,1);
    %%
%}

    %% Calculate the distance between each solution and point
%     V=0.4.*prod(max(PopDec)-min(PopDec)).^(1./size(PopDec,2));
    Distance = pdist2(PopDec,Wdec);
%     if Distance > V
%         Distance = 0;
%     end
%     Distance = mapminmax(Distance,0,1);
    Distanceobj = pdist2(PopObj,Wobj);%IGD
%     Distanceobj = mapminmax(Distanceobj,0,1);
    Con      = min(Distance,[],2);
    Cona      = min(Distanceobj,[],2);
    
%     Cona     = min(Distance,[],2);%IGD
    
%     %% 避免Con中出现全0
%     if sum(Con) == 0
%         for i = 1:N
%             AA = Temp(i, :);
%             AA(AA == 0) = inf;
%             Temp(i, :) = AA;
%         end
%         Con      = min(Temp,[],2);
%     end
%     %%
%     
    %% Delete the solution which has the smallest metric contribution one by one
    [dis,rank] = sort(Distance,1);
    [disa,ranka] = sort(Distanceobj,1);
    Remain     = true(1,N);
    while sum(Remain) > K
        % Calculate the fitness of outliers
        Outliers = Remain;
        t = Outliers;
        t(rank(1,:)) = false;
        tt = Outliers;
        tt(ranka(1,:)) = false;
        Outliers = tt & t;
        METRIC   = sum(dis(1,:)) + sum(Con(Outliers)) + sum(disa(1,:)) + sum(Cona(Outliers)); 
        Metrics  = inf(1,N);
        Metrics(Outliers) = METRIC - Con(Outliers)- Cona(Outliers);
        % Calculate the fitness of other solutions
        for p = find(Remain & ~Outliers)
            temp = rank(1,:) == p;
            temp1 = ranka(1,:) == p;
            outliers = false(1,N);
            outliers(rank(2,temp)) = true;
            outliers = outliers & Outliers;
            Metrics(p) = METRIC - sum(dis(1,temp)) + sum(dis(2,temp)) - sum(Con(outliers))- sum(disa(1,temp1)) + sum(disa(2,temp1)) - sum(Cona(outliers));
        end
        % Delete the worst solution and update the variables
        [~,del] = min(Metrics);
        temp = rank ~= del;
        temp1 = ranka ~= del;
        dis  = reshape(dis(temp),sum(Remain)-1,NW);
        rank = reshape(rank(temp),sum(Remain)-1,NW);
        disa  = reshape(disa(temp1),sum(Remain)-1,NW);
        ranka = reshape(ranka(temp1),sum(Remain)-1,NW);
        Remain(del) = false;
    end
end