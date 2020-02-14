function opt = TSEMO_options
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

% Description: Create TSEMO options structure.             
opt.maxeval = 35;                   % Maximum number of funciton evaluations

for i = 1:3                         % this is only valid for up to 3 objectives. Number needs to be increased for more than 3.
opt.GP(i).nSpectralpoints = 4000;   % Number of spectral sampling points for objective i
opt.GP(i).matern = 1;               % Matern type 1 / 3 / 5 / inf for objective i
opt.GP(i).fun_eval = 200;           % Function evaluations by direct algorithm per input dimension for objective i
end

opt.pop = 100;                      % Genetic algorithm population size
opt.Generation = 100;               % Genetic algorithm number of generations    
opt.NoOfBachSequential = 1;         % Batch size 
end
