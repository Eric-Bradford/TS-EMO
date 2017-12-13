function [pop, state] = evaluate(opt, pop, state, varargin)
% Function: [pop, state] = evaluate(opt, pop, state, varargin)
% Description: Evaluate the objective functions of each individual in the
%   population.
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************

N = length(pop);
allTime = zeros(N, 1);  % allTime : use to calculate average evaluation times

%*************************************************************************
% Evaluate objective function in parallel
%*************************************************************************
if( strcmpi(opt.useParallel, 'yes') == 1 )
    curPoolsize = matlabpool('size');

    % There isn't opened worker process
    if(curPoolsize == 0)
        if(opt.poolsize == 0)
            matlabpool open local
        else
            matlabpool(opt.poolsize)
        end
    % Close and recreate worker process
    else
        if(opt.poolsize ~= curPoolsize)
            matlabpool close
            matlabpool(opt.poolsize)
        end
    end

    parfor i = 1:N
        fprintf('\nEvaluating the objective function... Generation: %d / %d , Individual: %d / %d \n', state.currentGen, opt.maxGen, i, N);
        [pop(i), allTime(i)] = evalIndividual(pop(i), opt.objfun, varargin{:});
    end

%*************************************************************************
% Evaluate objective function in serial vectorized
%*************************************************************************
else
for i = 1:N
        fprintf('\nEvaluating the objective function... Generation: %d / %d , Individual: %d / %d \n', state.currentGen, opt.maxGen, i, N);
end
        [pop, allTime] = evalIndividual(pop, opt.objfun, varargin{:});
end

%*************************************************************************
% Statistics
%*************************************************************************
state.avgEvalTime   = sum(allTime) / length(allTime);
state.evaluateCount = state.evaluateCount + length(pop);




function [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% Function: [indi, evalTime] = evalIndividual(indi, objfun, varargin)
% Description: Evaluate one objective function.
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-25
%*************************************************************************

N = length(indi);
tStart = tic;
for i = 1:N
    pop_input(i,:) = indi(i).var;
end

[y, cons] = objfun(pop_input,varargin{:});
evalTime = toc(tStart);

% Save the objective values and constraint violations
for i = 1:N
indi(i).obj = y(i,:);
end



