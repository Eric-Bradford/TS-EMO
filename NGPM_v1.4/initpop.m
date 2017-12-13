function pop = initpop(opt, pop, varargin)
% Function: pop = initpop(opt, pop, varargin)
% Description: Initialize population.
% Syntax:
%   pop = initpop(opt, pop)
%     (default) Create a random initial population with a uniform distribution.
%
%   pop = initpop(opt, pop, 'pop.txt')
%     Load population from exist file and use the last population. If the popsize
%     less than the current popsize, then random numbers will used to fill the population. 
%
%   pop = initpop(opt, pop, 'pop.txt', ngen)
%     Load population from file with specified generation.
%
%   pop = initpop(opt, pop, oldresult)
%     Specify exist result structure.
%
%   pop = initpop(opt, pop, oldresult, ngen)
%     Specify exist result structure and the population which will be used.
%
% Parameters: 
%   pop : an empty population
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-01
%*************************************************************************


%*************************************************************************
% 1. Identify parameters
%*************************************************************************
method = 'uniform';

if(nargin >= 3)
    if( ischar(varargin{1}) )
        method = 'file';
    elseif( isstruct(varargin{1}) )
        method = 'existpop';
    end
end

%*************************************************************************
% 2. Initialize population with different methods
%*************************************************************************
if( strcmpi(method, 'uniform'))
    pop = initpopUniform(opt, pop);
elseif(strcmpi(method, 'file'))
    fprintf('...Initialize population from file "%s"\n', varargin{1});
    pop = initpopFromFile(opt, pop, varargin{:});
elseif(strcmpi(method, 'existpop'))
    fprintf('...Initialize population from specified result.\n');
    pop = initpopFromExistResult(opt, pop, varargin{:});
end




function pop = initpopFromFile(opt, pop, varargin)
% Function: pop = initpopFromFile(opt, pop, varargin)
% Description: Load population from specified population file.
% Syntax:
%   pop = initpop(opt, pop, 'pop.txt')
%   pop = initpop(opt, pop, 'pop.txt', ngen)
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************
fileName = varargin{1};

oldResult = loadpopfile(fileName);
pop = initpopFromExistResult(opt, pop, oldResult, varargin{2:end});





function pop = initpopFromExistResult(opt, pop, varargin)
% Function: pop = initpopFromExistResult(opt, pop, varargin)
% Description: Load population from exist result structure.
% Syntax:
%   pop = initpop(opt, pop, oldresult)
%   pop = initpop(opt, pop, oldresult, ngen)
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************

% 1. Verify param
oldresult = varargin{1};
if( ~isstruct(oldresult) || ~isfield(oldresult, 'pops') )
    error('NSGA2:InitPopError', 'The result structure specified is not correct!');
end


oldpops = oldresult.pops;
ind = oldpops(1,1);     % individual used to verify optimization param
if( opt.numVar ~= length(ind.var) || ...
    opt.numObj ~= length(ind.obj) || ...
    opt.numCons ~= length(ind.cons) )
    error('NSGA2:InitPopError', ...
        'The specified optimization result is not for current optimization model!');
end
clear ind


% 2. Determine which population would be used
ngen = 0;
if( nargin >= 4)
    ngen = varargin{2};
end

maxGen = size(oldpops, 1);
if(ngen == 0)
    ngen = maxGen;
elseif(ngen > maxGen)
    warning('NSGA2:InitPopWarning', ...
        'The specified generation "%d" does not exist, use "%d" instead.',...
        ngen, maxGen);
    ngen = maxGen;
end


% 3. Create initial population
popsizeOld = size(oldpops, 2);
popsizeNew = opt.popsize;

if( popsizeNew <= popsizeOld )      % a) All from old pop
    for i = 1:popsizeNew
        pop(i).var = oldpops(ngen, i).var;
    end
else                                % b) Use random individuals to fill the population
    for i = 1:popsizeOld
        pop(i).var = oldpops(ngen, i).var;
    end
    pop(popsizeOld+1:end) = initpopUniform(opt, pop(popsizeOld+1:end));
end




function pop = initpopUniform(opt, pop)
% Function: pop = initpopUniform(opt, pop)
% Description: Initialize population using random number
%
%    Copyright 2011 by LSSSSWC
%    Revision: 1.0  Data: 2011-07-01
%*************************************************************************

nVar = opt.numVar;
type = opt.vartype;

lb = opt.lb;
ub = opt.ub;

popsize = length(pop);
for i = 1:popsize
    var = lb + rand(1, nVar) .* (ub-lb);
    
    % if desing variable is integer, round to the nearest integer
    for v = 1:nVar
        if( type(v) == 2)
            var(v) = round(var(v));
        end
    end
    
    % limit in the lower and upper bound
    var = varlimit(var, lb, ub);

    pop(i).var = var;

end






