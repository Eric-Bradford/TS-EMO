function pop = mutationOp(opt, pop, state)
% Function: pop = mutationOp(opt, pop, state)
% Description: Mutation Operator. All of the individuals would do mutation, but
%   only "mutationFraction" of design variables of an individual would changed.
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************

%*************************************************************************
% 1. Check for the parameters
%*************************************************************************
% mutation method
strfun = lower(opt.mutation{1});
numOptions = length(opt.mutation) - 1;
[mutationopt{1:numOptions}] = opt.mutation{2:end};

switch (strfun)
    case 'gaussian'
        fun = @mutationGaussian;
    otherwise
        error('NSGA2:MutationOpError', 'No support mutation operator!');
end

nVar = opt.numVar;

% "auto" mutation fraction
if( ischar(opt.mutationFraction) )
    if( strcmpi(opt.mutationFraction, 'auto') )
        fraction = 2.0 / nVar;
    else
        error('NSGA2:MutationOpError', 'The "mutationsFraction" parameter should be scalar or "auto" string.');
    end
else
    fraction = opt.mutationFraction;
end


% All of the individual would be modified, but only 'mutationFraction' of design
% variables for an individual would be changed.
for ind = 1:length(pop)
        child = fun( pop(ind), opt, state, fraction, mutationopt);
        
        % Rounding for integer variables
        for v = 1:nVar
            if( opt.vartype(v) == 2)
                child.var(v) = round( child.var(v) );
            end
        end

        child.var = varlimit(child.var, opt.lb, opt.ub);
        
        pop(ind) = child;
end



function child = mutationGaussian( parent, opt, state, fraction, options)
% Function: child = mutationGaussian( parent, opt, state, fraction, options)
% Description: Gaussian mutation operator. Reference Matlab's help :
%   Genetic Algorithm Options :: Options Reference (Global Optimization Toolbox)
% Parameters: 
%   fraction : mutation fraction of variables of an individual
%   options{1} : scale. This paramter should be large enough for interger variables
%     to change from one to another.
%   options{2} : shrink
% Return: 
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-13
%*************************************************************************


%*************************************************************************
% 1. Verify the parameters.
%*************************************************************************
if( length(options)~=2)
    error('NSGA2:MutationOpError', 'Mutation operator parameter error!');
end


%*************************************************************************
% 2. Calc the "scale" and "shrink" parameter.
%*************************************************************************
scale = options{1};
shrink = options{2};
scale = scale - shrink * scale * state.currentGen / opt.maxGen;

lb = opt.lb;
ub = opt.ub;
scale = scale * (ub - lb);


%*************************************************************************
% 3. Do the mutation.
%*************************************************************************
child = parent;
numVar = length(child.var);
for i = 1:numVar
    if(rand() < fraction)
        child.var(i) = parent.var(i) + scale(i) * randn();
    end
end





