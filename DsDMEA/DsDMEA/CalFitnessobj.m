function Fitness = CalFitnessobj(PopObj,PopDec)
% Calculate the fitness of each solution

N = size(PopObj,1);

%% Detect the dominance relation between each two solutions
Dominate = false(N);
for i = 1 : N-1
    for j = i+1 : N
        k = any(PopObj(i,:)<PopObj(j,:)) - any(PopObj(i,:)>PopObj(j,:));
        if k == 1
            Dominate(i,j) = true;
        elseif k == -1
            Dominate(j,i) = true;
        end
    end
end


dist=pdist2(PopDec,PopDec);
V=0.5.*prod(max(PopDec)-min(PopDec)).^(1./size(PopDec,2));
S = zeros(N,1);
for i=1:N
    index = (dist(i,:)<V);
    S(i) = sum(Dominate(i,index));
end


K = zeros(1,N);
for i = 1 : N
    index = (dist(i,:)<V);
    K(i) = sum(S(Dominate(index,i)));
end


Distance = pdist2(PopObj,PopObj);
Distance(logical(eye(length(Distance)))) = inf;
Distance = sort(Distance,2);
D = 1./(Distance(:,floor(sqrt(N)))+2);

Fitness = K + D';
end