function [opt, pop] = ndsort(opt, pop)
% Function: [opt, pop] = ndsort(pop)
% Description: Fast non-dominated sort.
%
%         LSSSSWC, NWPU
%    Revision: 1.4  Data: 2011-07-26
%*************************************************************************




%*************************************************************************
% 1. Initialize variables
%   indi.np£ºnumber of individuals which dominate this individual
%   indi.sp(:): a set of individuals that this individual dominate
%*************************************************************************
N = length(pop);    %popsize
ind = repmat(struct('np',0, 'sp', []),[1,N]);

for i = 1:N
    pop(i).rank = 0;
    pop(i).distance = 0;
    pop(i).prefDistance = 0;
end


%*************************************************************************
% 2. fast non-dominated sort
%*************************************************************************
% Calculate the domination matrix for improving the efficiency.

% NOTE: The "for" statement is more efficient than "vertcat" statement in my computer
% on Matlab 2010b. I don't know why.(LSSSSWC, 2011-07-25)
nViol   = zeros(N, 1);
violSum = zeros(N, 1);
for i = 1:N
    nViol(i)    = pop(i).nViol;
    violSum(i)  = pop(i).violSum;
end
% nViol   = vertcat(pop(:).nViol);
% violSum = vertcat(pop(:).violSum);

obj     = vertcat(pop(:).obj);
domMat  = calcDominationMatrix(nViol, violSum, obj); % domination matrix for efficiency


% Compute np and sp of each indivudal
for p = 1:N-1
    for q = p+1:N
        if(domMat(p, q) == 1)          % p dominate q
            ind(q).np = ind(q).np + 1;
            ind(p).sp = [ind(p).sp , q];
        elseif(domMat(p, q) == -1)     % q dominate p
            ind(p).np = ind(p).np + 1;
            ind(q).sp = [ind(q).sp , p];
        end
    end
end


% The first front(rank = 1)
front(1).f = [];    % There are only one field 'f' in structure 'front'.
                    % This is intentional because the number of individuals
                    % in the front is difference.
for i = 1:N
    if( ind(i).np == 0 )
        pop(i).rank = 1;
        front(1).f = [front(1).f, i];
    end
end

% Calculate pareto rank of each individuals, viz., pop(:).rank 
fid = 1;        %pareto front ID
while( ~isempty(front(fid).f) )
    Q = [];
    for p = front(fid).f
        for q = ind(p).sp
            ind(q).np = ind(q).np -1;
            if( ind(q).np == 0 )
                pop(q).rank = fid+1;
                Q = [Q, q];
            end
        end
    end
    fid = fid + 1;
    
    front(fid).f = Q;
end
front(fid) = [];    % delete the last empty front set



%*************************************************************************
% 3. Calculate the distance
%*************************************************************************
if(isempty(opt.refPoints))
    pop = calcCrowdingDistance(opt, pop, front);
else
    [opt, pop] = calcPreferenceDistance(opt, pop, front);
end





function domMat = calcDominationMatrix(nViol, violSum, obj)
% Function: domMat = calcDominationMatrix(nViol, violSum, obj)
% Description: Calculate the domination maxtir which specified the domination
%   releation between two individual using constrained-domination.
%
% Return: 
%   domMat(N,N) : domination matrix
%       domMat(p,q)=1  : p dominates q
%       domMat(p,q)=-1 : q dominates p
%       domMat(p,q)=0  : non dominate
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-13
%*************************************************************************

N       = size(obj, 1);
numObj  = size(obj, 2);

domMat  = zeros(N, N);

for p = 1:N-1
    for q = p+1:N
        %*************************************************************************
        % 1. p and q are both feasible
        %*************************************************************************
        if(nViol(p) == 0 && nViol(q)==0)
            pdomq = false;
            qdomp = false;
            for i = 1:numObj
                if( obj(p, i) < obj(q, i) )         % objective function is minimization!
                    pdomq = true;
                elseif(obj(p, i) > obj(q, i))
                    qdomp = true;
                end
            end

            if( pdomq && ~qdomp )
                domMat(p, q) = 1;
            elseif(~pdomq && qdomp )
                domMat(p, q) = -1;
            end
        %*************************************************************************
        % 2. p is feasible, and q is infeasible
        %*************************************************************************
        elseif(nViol(p) == 0 && nViol(q)~=0)
            domMat(p, q) = 1;
        %*************************************************************************
        % 3. q is feasible, and p is infeasible
        %*************************************************************************
        elseif(nViol(p) ~= 0 && nViol(q)==0)
            domMat(p, q) = -1;
        %*************************************************************************
        % 4. p and q are both infeasible
        %*************************************************************************
        else
            if(violSum(p) < violSum(q))
                domMat(p, q) = 1;
            elseif(violSum(p) > violSum(q))
                domMat(p, q) = -1;
            end
        end
    end
end

domMat = domMat - domMat';





function [opt, pop] = calcPreferenceDistance(opt, pop, front)
% Function: [opt, pop] = calcPreferenceDistance(opt, pop, front)
% Description: Calculate the 'preference distance' used in R-NSGA-II.
% Return: 
%   opt : This structure may be modified only when opt.refUseNormDistance=='ever'.
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.1  Data: 2011-07-26
%*************************************************************************

%*************************************************************************
% 1. Initialization
%*************************************************************************
numObj = length( pop(1).obj );  % number of objectives

refPoints = opt.refPoints;
refWeight = opt.refWeight;      % weight factor of objectives
if(isempty(refWeight))
    refWeight = ones(1, numObj);
end
epsilon = opt.refEpsilon;
numRefPoint = size(refPoints, 1);

% Determine the normalized factors
bUseFrontMaxMin = false;    % bUseFrontMaxMin : If use the maximum and minimum value in the front as normalized factor.
if( strcmpi(opt.refUseNormDistance, 'ever') )
    % 1) Find possiable (not current population) maximum and minimum value
    %     of each objective.
    obj = vertcat(pop.obj);
    if( ~isfield(opt, 'refObjMax_tmp') )
        opt.refObjMax_tmp = max(obj);
        opt.refObjMin_tmp = min(obj);
    else
        objMax = max(obj);
        objMin = min(obj);
        for i = 1:numObj
            if(opt.refObjMax_tmp(i) < objMax(i))
                opt.refObjMax_tmp(i) = objMax(i);
            end
            if(opt.refObjMin_tmp(i) > objMin(i))
                opt.refObjMin_tmp(i) = objMin(i);
            end
        end
        clear objMax objMin
    end
    objMaxMin = opt.refObjMax_tmp - opt.refObjMin_tmp;
    clear obj
elseif( strcmpi(opt.refUseNormDistance, 'front') )
    % 2) Do not use normalized Euclidean distance.
    bUseFrontMaxMin = true;
elseif( strcmpi(opt.refUseNormDistance, 'no') )
    % 3) Do not use normalized Euclidean distance.
    objMaxMin = ones(1,numObj);
else
    % 3) Error
    error('NSGA2:ParamError', ...
        'No support parameter: options.refUseNormDistance="%s", only "yes" or "no" are supported',...
        opt.refUseNormDistance);
end


%*************************************************************************
% 2. Calculate preference distance pop(:).prefDistance
%*************************************************************************
for fid = 1:length(front)
    % Step1: Calculate the weighted Euclidean distance in each front
    idxFront = front(fid).f;            % idxFront : index of individuals in current front
    numInd = length(idxFront);          % numInd : number of individuals in current front
    popFront = pop(idxFront);           % popFront : individuals in front fid

    objFront = vertcat(popFront.obj);   % objFront : the whole objectives of all individuals

    if(bUseFrontMaxMin)
        objMaxMin = max(objFront) - min(objFront); % objMaxMin : the normalized factor in current front
    end

    % normDistance : weighted normalized Euclidean distance
    normDistance = calcWeightNormDistance(objFront, refPoints, objMaxMin, refWeight);
    
    
    % Step2: Assigned preference distance
    prefDistanceMat = zeros(numInd, numRefPoint);
    for ipt = 1:numRefPoint
        [~,ix] = sort(normDistance(:, ipt));
        prefDistanceMat(ix, ipt) = 1:numInd;
    end
    prefDistance = min(prefDistanceMat, [], 2);
    clear ix

    
    % Step3: Epsilon clearing strategy
    idxRemain = 1:numInd;           % idxRemain : index of individuals which were not processed
    while(~isempty(idxRemain))
        % 1. Select one individual from remains
        objRemain = objFront( idxRemain, :);
        selIdx = randi( [1,length(idxRemain)] );
        selObj = objRemain(selIdx, :);

        % 2. Calc normalized Euclidean distance
        % distanceToSel : normalized Euclidean distance to the selected points
        distanceToSel = calcWeightNormDistance(objRemain, selObj, objMaxMin, refWeight);
        

        % 3. Process the individuals within a epsilon-neighborhood
        idx = find( distanceToSel <= epsilon );     % idx : index in idxRemain
        if(length(idx) == 1)    % the only individual is the selected one
            idxRemain(selIdx)=[];
        else
            for i=1:length(idx)
                if( idx(i)~=selIdx )
                    idInIdxRemain = idx(i);     % idx is the index in idxRemain vector
                    id = idxRemain(idInIdxRemain);
                    
                    % *Increase the preference distance to discourage the individuals
                    % to remain in the selection.
                    prefDistance(id) = prefDistance(id) + round(numInd/2);
                end
            end
            idxRemain(idx) = [];
        end
        
    end

    % Save the preference distance
    for i=1:numInd
        id = idxFront(i);
        pop(id).prefDistance = prefDistance(i);
    end
end


function distance = calcWeightNormDistance(points, refPoints, maxMin, weight)
% Function: calcWeightNormDistance(points, refPoints, maxMin, weight)
% Description: Calculate the weighted Euclidean distance from "points" to "refPoints"
% Parameters: 
%   points(nPoint, N)       : each row is a point in N dimension space.
%   refPoints(nRefPoint, N) : each row is a reference point.
%   maxMin(1, N)            : normalized factor.
%   weight(1, N)            : weights
%
% Return: 
%   distance(nPoint, nRefPoint)
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-14
%*************************************************************************

nRefPoint = size(refPoints, 1);     % number of reference points
nPoint = size(points, 1);           % number of points

distance = zeros(nPoint, nRefPoint);
for ipt = 1:nRefPoint
    refpt = refPoints(ipt, :);
    for i = 1:nPoint
        weightNormDist = ((points(i, :)-refpt) ./ maxMin).^2 .* weight;
        distance(i, ipt) = sqrt(sum(weightNormDist));
    end
end





function pop = calcCrowdingDistance(opt, pop, front)
% Function: pop = calcCrowdingDistance(opt, pop, front)
% Description: Calculate the 'crowding distance' used in the original NSGA-II.
% Syntax:
% Parameters: 
% Return: 
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-11
%*************************************************************************

numObj = length( pop(1).obj );  % number of objectives
for fid = 1:length(front)
    idx = front(fid).f;
    frontPop = pop(idx);        % frontPop : individuals in front fid
    
    numInd = length(idx);       % nInd : number of individuals in current front
    
    obj = vertcat(frontPop.obj);
    obj = [obj, idx'];          % objctive values are sorted with individual ID
    for m = 1:numObj
        obj = sortrows(obj, m);

        colIdx = numObj+1;
        pop( obj(1, colIdx) ).distance = Inf;         % the first one
        pop( obj(numInd, colIdx) ).distance = Inf;    % the last one
        
        minobj = obj(1, m);         % the maximum of objective m
        maxobj = obj(numInd, m);    % the minimum of objective m
        
        for i = 2:(numInd-1)
            id = obj(i, colIdx);
            pop(id).distance = pop(id).distance + (obj(i+1, m) - obj(i-1, m)) / (maxobj - minobj);
        end
    end
end






