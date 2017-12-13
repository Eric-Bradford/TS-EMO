function pop = crossoverOp(opt, pop, state)
% Function: pop = crossoverOp(opt, pop, state)
% Description: Crossover operator. All of the individuals would be do crossover, but
%   only "crossoverFraction" of design variables of an individual would changed.
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************

%*************************************************************************
% 1. Check for the parameters
%*************************************************************************
% determine the crossover method
strfun = lower(opt.crossover{1});
numOptions = length(opt.crossover) - 1;
[crossoverOpt{1:numOptions}] = opt.crossover{2:end};

switch( strfun )
    case 'intermediate'
        fun = @crsIntermediate;
    otherwise
        error('NSGA2:CrossoverOpError', 'No support crossover operator!');
end

nVar = opt.numVar;

% "auto" crossover fraction
if( ischar(opt.crossoverFraction) )
    if( strcmpi(opt.crossoverFraction, 'auto') )
        fraction = 2.0 / nVar;
    else
        error('NSGA2:CrossoverOpError', 'The "crossoverFraction" parameter should be scalar or "auto" string.');
    end
else
    fraction = opt.crossoverFraction;
end


for ind = 1:2:length(pop)    % Popsize should be even number
        % Create children
        [child1, child2] = fun( pop(ind), pop(ind+1), fraction, crossoverOpt );
        
        % Round
        for v = 1:nVar
            if( opt.vartype(v) == 2)
                child1.var(v) = round( child1.var(v) );
                child2.var(v) = round( child2.var(v) );
            end
        end

        % Bounding limit
        child1.var = varlimit(child1.var, opt.lb, opt.ub);
        child2.var = varlimit(child2.var, opt.lb, opt.ub);
        
        pop(ind)     = child1;
        pop(ind+1)   = child2;
    
end



function [child1, child2] = crsIntermediate(parent1, parent2, fraction, options)
% Function: [child1, child2] = crsIntermediate(parent1, parent2, fraction, options)
% Description: (For real coding) Intermediate crossover. (Same as Matlab's crossover 
%   operator)
%       child = parent1 + rand * Ratio * ( parent2 - parent1)
% Parameters: 
%   fraction : crossover fraction of variables of an individual
%   options = ratio
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************


if( length(options)~=1 || ~isnumeric(options{1}))
    error('NSGA2:CrossoverOpError', 'Crossover operator parameter error!');
end

ratio = options{1};

child1 = parent1;
child2 = parent2;

nVar = length(parent1.var);
crsFlag = rand(1, nVar) < fraction;

randNum = rand(1,nVar);     % uniformly distribution

child1.var = parent1.var + crsFlag .* randNum .* ratio .* (parent2.var - parent1.var);
child2.var = parent2.var - crsFlag .* randNum .* ratio .* (parent2.var - parent1.var);




