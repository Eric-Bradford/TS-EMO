function [Xpareto,Ypareto,X,Y,XParetoGP,YParetoGP,hypf] = TSEMO_V3(f,X,Y,lb,ub,Opt)
%  Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2018-08-07.
%  Version 3 by Eric Bradford (eric.bradford@ntnu.no), Artur M.
%  Schweidtmann (Artur.Schweidtmann@avt.rwth-aachen.de) and Alexei Lapkin
%  (aal35@cam.ac.uk)

% ALGORITHM  Executes the optimization algorithm.

%   OUTPUTS
%       Xpareto     Pareto front approximation from X inputs       [np,D]
%       Ypareto     Pareto front approximation from Y outputs      [np,O]
%       X           Data input vector                              [n+meval,D]
%       Y           Data output vector                             [n+meval,O]
%       XParetoGP   Pareto front approximation of GP model inputs  [npGP,D]
%       YParetoGP   Pareto front approximation of GP model outputs [npGP,D]
%       hypf        Final hyperparameter values []                 [D+2,O]

%   INPUTS
%       f       objective function
%       X       input values of given data set  [n,D]
%       Y       output values of given data set [n,O]
%       lb      lower bounds of input variables [1,D]
%       ub      upper bounds of input variables [1,D]
%       Opt     options struct of algorithm

%   with
%       n       number of given data points
%       meval   maximum number of function evaluations
%       np      number of pareto points of the final dataset
%       npGP    number of pareto points of GP model dataset
%       D       dimension of input space
%       O       dimension of output space

%% Initialize option structure
it = 1;
Opt = set_option_structure(Opt,X,Y);

%% Write initial text file for TSEMO_log
create_log_file(X,Y,Opt,f,lb,ub)

for i = 1:ceil(Opt.maxeval/Opt.NoOfBachSequential)
    tic;
    %% Scale Variables
    [Xnew,Ynew] = ScaleVariables(X,Y,lb,ub,Opt) ;
    
    %% Training of GP
    for j = 1:Opt.Gen.NoOfGPs
        Opt.GP(j).hyp = TrainingOfGP(Xnew,Ynew(:,j),Opt.GP(j));
    end
    
    %% Draw samples of GPs
    for j = 1:Opt.Gen.NoOfGPs
        Opt.Sample(j).f = posterior_sample(Xnew,Ynew(:,j),Opt.GP(j));
    end
    
    %% Determine pareto front of function samples
    [Sample_pareto,Sample_xpareto,Sample_nadir] = Find_sample_pareto(Opt);
    
    %% Sampling point at maximum hypervolume improvement
    [index,hv_imp] = hypervolume_improvement_index(Ynew,Sample_nadir,Sample_pareto,Opt);
    
    xNew = Sample_xpareto(index,:);
    Xnew = [Xnew;xNew];
    
    for j = 1:Opt.Gen.NoOfInputDim
        xnewtrue(:,j) = xNew(:,j)*(ub(j)-lb(j)) + lb(j);
    end
    
    for l = 1 : size(xnewtrue,1)
        ytrue(l,:) = f(xnewtrue(l,:));
    end
    
    X = [X;xnewtrue];
    Y = [Y;ytrue];
    Ynew = zeros(size(Y,1),Opt.Gen.NoOfGPs);
    for j = 1:Opt.Gen.NoOfGPs
        Ynew(:,j) = (Y(:,j) - mean(Y(:,j)))/std(Y(:,j));
    end
    
    front = paretofront(Y);
    Xpareto = X(front,:);
    Ypareto = Y(front,:);
    
    %% Update log each iterations
    update_log_file(it,hv_imp,toc,xnewtrue,ytrue,Opt,Y,ub,lb)
    
    %% Determine final hyperparameter values for analysis
    if i == ceil(Opt.maxeval/Opt.NoOfBachSequential)
        for j = 1:Opt.Gen.NoOfGPs
            Opt.GP(j).hyp = TrainingOfGP(Xnew,Ynew(:,j),Opt.GP(j));
        end
        
        hypf = zeros(Opt.Gen.NoOfInputDim+2,Opt.Gen.NoOfGPs);
        for j = 1:Opt.Gen.NoOfInputDim
            covhyp = exp(Opt.GP(j).hyp.cov);
            hypf(:,j) = [covhyp(1:Opt.Gen.NoOfInputDim)*(1/(ub(j)-lb(j)));covhyp(end)*std(Y(:,j));exp(Opt.GP(j).hyp.lik)*std(Y(:,j))];
        end
        
        %% Obtain Pareto front from spectral Gaussian process model
        for j = 1:Opt.Gen.NoOfGPs
            Opt.Mean(j).f = mean_sample(Xnew,Ynew(:,j),Opt.GP(j));
        end
        
        [Mean_pareto,Mean_xpareto] = Find_mean_pareto(Opt);
        
        for j = 1:Opt.Gen.NoOfInputDim
            XParetoGP(:,j) = Mean_xpareto(:,j)*(ub(j)-lb(j)) + lb(j);
        end
        
        for j = 1:Opt.Gen.NoOfGPs
            YParetoGP(:,j) = Mean_pareto(:,j)*std(Y(:,j)) + mean(Y(:,j));
        end
        
        
        %% Update log with final results
        final_log_update(Xpareto,Ypareto,X,Y,XParetoGP,YParetoGP,hypf,Opt)
    end
    
    %% Display
    if it == 1
        fprintf('%10s %10s %10s \n','Iteration', 'HypImp', 'Time(s)');
    end
    fprintf('%10d %10.4g %10.3g \n', it, hv_imp, toc);
    
    it = it+1;
    
end
return

function create_log_file(X,Y,Opt,f,lb,ub)
try
    function_name = func2str(f);
    string1 = '';
    string2 = {};
    string3 = '';
    string4 = '';
    string5 = {};
    string6 = '';
    string7 = '';
    for i = 1:size(X,2)
        string1 = strcat(string1,'%8.4f ');
        string2 = {string2{:},strcat('x',num2str(i))};
        string3 = strcat(string3,'%+8s ');
    end
    for i = 1:size(Y,2)
        string4 = strcat(string4,'%8.4f ');
        string5 = {string5{:},strcat('f',num2str(i))};
        string6 = strcat(string6,'%+8s ');
        string7 = strcat(string7,'%8d ');
    end
    
    TSEMO_log = fopen( 'TSEMO_log.txt', 'w');
    fprintf(TSEMO_log,'\n %s %s \n', 'TSEMO log file created on',date);
    fprintf(TSEMO_log,'\n %s','This file shows the initial specifications of TSEMO and logs the output.');
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'License information');
    fprintf(TSEMO_log,'\n %s \n', 'BSD 2-Clause License');
    fprintf(TSEMO_log,'\n %s', 'Copyright (c) 2017, Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin');
    fprintf(TSEMO_log,'\n %s \n', 'All rights reserved.');
    fprintf(TSEMO_log,'\n %s', 'Redistribution and use in source and binary forms, with or without');
    fprintf(TSEMO_log,'\n %s \n', 'modification, are permitted provided that the following conditions are met:');
    fprintf(TSEMO_log,'\n %s   ', '*Redistributions of source code must retain the above copyright notice, this');
    fprintf(TSEMO_log,'\n %s \n', ' list of conditions and the following disclaimer.');
    fprintf(TSEMO_log,'\n %s   ', '*Redistributions in binary form must reproduce the above copyright notice,');
    fprintf(TSEMO_log,'\n %s   ', ' this list of conditions and the following disclaimer in the documentation');
    fprintf(TSEMO_log,'\n %s \n', ' and/or other materials provided with the distribution.');
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'Problem specifications');
    fprintf(TSEMO_log,'\n %s %s \n', 'Function used:  ',function_name);
    fprintf(TSEMO_log,'\n %s %d', 'Number of inputs:  ',size(X,2));
    fprintf(TSEMO_log,'\n %s %d \n', 'Number of outputs: ',size(Y,2));
    fprintf(TSEMO_log,'\n %s', 'Lower bounds of decision variables:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,string1,lb);
    fprintf(TSEMO_log,'\n \n %s', 'Upper bounds of decision variables:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,string1,ub);
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'Algorithm options');
    fprintf(TSEMO_log,'\n %s %d', 'Maximum number of function evaluations: ',Opt.maxeval);
    fprintf(TSEMO_log,'\n %s %d', 'Sample batch size:                      ',Opt.NoOfBachSequential);
    fprintf(TSEMO_log,'\n %s %d \n', 'Number of algorithm iterations:         ',ceil(Opt.maxeval/Opt.NoOfBachSequential));
    fprintf(TSEMO_log,'\n %s %d', 'Genetic algorithm population size:       ',Opt.pop);
    fprintf(TSEMO_log,'\n %s %d \n', 'Genetic algorithm number of generations: ',Opt.Generation);
    fprintf(TSEMO_log,strcat('\n','%s',string6,'\n'),'                                         ',string5{:});
    fprintf(TSEMO_log,strcat('%s',string7,'\n'), ' Number of spectral sampling points:     ',Opt.GP(1:size(Y,2)).nSpectralpoints);
    fprintf(TSEMO_log,strcat('%s',string7), ' Type of matern function:                ',Opt.GP(1:size(Y,2)).matern);
    fprintf(TSEMO_log,strcat('\n','%s',string7), ' Direct evaluations per input dimension: ',Opt.GP(1:size(Y,2)).fun_eval);
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    
    fprintf(TSEMO_log,'\n %s \n', 'Initial data set');
    fprintf(TSEMO_log,'\n %s %d \n', 'Number of initial data points: ',size(X,1));
    fprintf(TSEMO_log,'\n %s', 'Initial input data matrix:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),X');
    fprintf(TSEMO_log,'\n %s', 'Initial output data matrix:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),Y');
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    fclose(TSEMO_log) ;
catch
    warning('There was an error when writing the log file. Maybe the file was currently open in another program. The algorithm continues to run but some parts of the log file are maybe not correct.') ;
    fclose('all') ;
end
return

function update_log_file(it,hv_imp,toc,xnewtrue,ytrue,Opt,Y,ub,lb)
try
    string1 = '';
    string2 = {};
    string3 = '';
    string4 = '';
    string5 = {};
    string6 = '';
    string7 = '';
    string8 = {};
    for i = 1:Opt.Gen.NoOfInputDim
        string1 = strcat(string1,'%8.4f ');
        string2 = {string2{:},strcat('x',num2str(i))};
        string3 = strcat(string3,'%+8s ');
        string8 = {string8{:},strcat('lambda',num2str(i))};
    end
    string8 = {string8{:},'sigmaf'};
    string8 = {string8{:},'sigman'};
    for i = 1:Opt.Gen.NoOfGPs
        string4     = strcat(string4,'%8.4f ');
        string5     = {string5{:},strcat('f',num2str(i))};
        string6     = strcat(string6,'%+8s ');
        string7     = strcat(string7,'%8d ');
        hypcov      = exp(Opt.GP(i).hyp.cov);
        hypmat(:,i) = [hypcov(1:end-1).*(1./(ub-lb))';hypcov(end)*std(Y(:,i));exp(Opt.GP(i).hyp.lik)*std(Y(:,i))];
    end
    
    TSEMO_log = fopen( 'TSEMO_log.txt', 'a');
    fprintf(TSEMO_log,'\n %s %d \n', 'Algorithm iteration',it);
    fprintf(TSEMO_log,'\n %s %8.4f', 'Predicted hypervolume improvement: ',hv_imp);
    fprintf(TSEMO_log,'\n %s %8.4f \n', 'Time taken: ',toc);
    fprintf(TSEMO_log,'\n %s', 'Proposed evaluation point(s): ');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),xnewtrue');
    fprintf(TSEMO_log,'\n %s', 'Corresponding observation(s): ');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),ytrue);
    fprintf(TSEMO_log,'\n %s', 'Current hyperparameter values: ');
    fprintf(TSEMO_log,strcat('\n','%+16s',string6,'\n'),'Hyperparameter',string5{:});
    for i = 1:Opt.Gen.NoOfInputDim+2
        fprintf(TSEMO_log,strcat('%+16s',string4,'\n'),string8{i},hypmat(i,:));
    end
    fprintf(TSEMO_log,'\n %s', '¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯');
    fclose(TSEMO_log) ;
catch
    warning('There was an error when writing the log file. Maybe the file was currently open in another program. The algorithm continues to run but some parts of the log file are maybe not correct.') ;
    fclose('all') ;
end
return

function final_log_update(Xpareto,Ypareto,X,Y,XParetoGP,YParetoGP,hypf,Opt)
try
    string1 = '';
    string2 = {};
    string3 = '';
    string4 = '';
    string5 = {};
    string6 = '';
    string7 = '';
    string8 = {};
    for i = 1:Opt.Gen.NoOfInputDim
        string1 = strcat(string1,'%8.4f ');
        string2 = {string2{:},strcat('x',num2str(i))};
        string3 = strcat(string3,'%+8s ');
        string8 = {string8{:},strcat('lambda',num2str(i))};
    end
    string8 = {string8{:},'sigmaf'};
    string8 = {string8{:},'sigman'};
    for i = 1:Opt.Gen.NoOfGPs
        string4     = strcat(string4,'%8.4f ');
        string5     = {string5{:},strcat('f',num2str(i))};
        string6     = strcat(string6,'%+8s ');
        string7     = strcat(string7,'%8d ');
    end
    
    TSEMO_log = fopen( 'TSEMO_log.txt', 'a');
    fprintf(TSEMO_log,'\n %s \n', 'Final algorithm output');
    fprintf(TSEMO_log,'\n %s', 'Final input data matrix:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),X');
    fprintf(TSEMO_log,'\n %s', 'Final output data matrix:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),Y');
    fprintf(TSEMO_log,'\n %s', 'Input data matrix of corresponding Pareto front:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),Xpareto');
    fprintf(TSEMO_log,'\n %s', 'Output data matrix of corresponding Pareto front:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),Ypareto');
    fprintf(TSEMO_log,'\n %s', 'Input data matrix of Pareto front of final Gaussian process model:');
    fprintf(TSEMO_log,strcat('\n',string3,'\n'),string2{:});
    fprintf(TSEMO_log,strcat(string1,'\n'),XParetoGP');
    fprintf(TSEMO_log,'\n %s', 'Output data matrix of Pareto front of final Gaussian process model:');
    fprintf(TSEMO_log,strcat('\n',string6,'\n'),string5{:});
    fprintf(TSEMO_log,strcat(string4,'\n'),YParetoGP');
    fprintf(TSEMO_log,'\n %s', 'Final hyperparameter values: ');
    fprintf(TSEMO_log,strcat('\n','%+16s',string6,'\n'),'Hyperparameter',string5{:});
    for i = 1:Opt.Gen.NoOfInputDim+2
        fprintf(TSEMO_log,strcat('%+16s',string4,'\n'),string8{i},hypf(i,:));
    end
    fclose(TSEMO_log) ;
catch
    warning('There was an error when writing the log file. Maybe the file was currently open in another program. The algorithm continues to run but some parts of the log file are maybe not correct.') ;
    fclose('all') ;
end
return

function [Xnew,Ynew] = ScaleVariables(X,Y,lb,ub,Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

%% Scales input and output variabels
Xnew = zeros(size(X)) ; % scaled inputs
Ynew = zeros(size(Y)) ; % scaled outputs

%% Scale input variables to [0,1]
for i = 1 : size(X,2)
    Xnew(:,i) = (X(:,i)-lb(i)) / (ub(i)-lb(i)) ;
end

%% Scale output variables to zero mean and unit variance
MeanOfOutputs = zeros(Opt.Gen.NoOfGPs,1) ;
stdOfOutputs = zeros(Opt.Gen.NoOfGPs,1) ;
for i = 1 : size(Y,2)
    MeanOfOutputs(i) = mean(Y(:,i)); % calculate mean
    stdOfOutputs(i) = std(Y(:,i)) ; % calculate standard deviation
    Ynew(:,i) = (Y(:,i) - MeanOfOutputs(i)) / stdOfOutputs(i) ; % scale outputs
end
return

function [OptGPhyp] = TrainingOfGP(Xnew,Ynew,OptGP)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

% Function which minimizes the neg-loglikelihood to find hyperparameters
%% Initialize variables
Opt.GP = OptGP ;
[n,D] = size(Xnew) ;

%% Set initial hyperparameters
h1 = Opt.GP.h1; % number of hyperparameters from covariance
h2 = Opt.GP.h2; % number of hyperparameters from likelihood

%% Calculation of squared-distance matrix
a = Xnew' ;
K_M = zeros(n,n*D) ;
for i = 1:D
    K_M(:,(i-1)*n+1:i*n) = sqdist(a(i,:),a(i,:)) ;
end

%% Minimize log-negative likeliehood
% Objective Function
obj_fun.f = @(hypVar) NLikelihood(hypVar,Xnew,Ynew,K_M,Opt.GP);

% Define bounds
lb              = ones(h1+h2,1) *  log(sqrt(10^(-3))) ;  % see Jones paper
ub              = ones(h1+h2,1) *  log(sqrt(10^(3)))  ;  % see Jones paper
lb(h1+h2)    = -6;
ub(h1+h2)    = Opt.GP.noiselimit;
bounds = [lb,ub];
opts.maxevals = Opt.GP.fun_eval*(h1+h2);
opts.maxits =  100000*(h1+h2);
opts.maxdeep = 100000*(h1+h2);
opts.showits = 0;

% Defintion of options for global search
[~,x0] = Direct(obj_fun,bounds,opts);

% Defintion of options for fmincon solver
LSoptions.Algorithm = 'interior-point';
LSoptions.DerivativeCheck = 'off';
LSoptions.TolCon = 1e-12;
LSoptions.Display = 'off';
LSoptions.Hessian = 'bfgs';
LSoptions.TolFun = 1e-12;
LSoptions.PlotFcns = [];
LSoptions.GradConstr = 'off';
LSoptions.GradObj = 'on';
LSoptions.TolX = 1e-14;
LSoptions.UseParallel = 0;

% Solve optimization problem
hypResult = fmincon(obj_fun.f,x0,[],[],[],[],lb,ub,[],LSoptions);

%% Return optimal hyperparameters
Opt.GP.hyp.cov  = hypResult(1:h1);
Opt.GP.hyp.lik  = hypResult(h1+1:h1+h2);
OptGPhyp = Opt.GP.hyp ;
return

function [NLL,dNLL] = NLikelihood(hypVar, Xnew, Ynew, K_M, OptGP)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

% Calculates the log-negative likelihood
%% Initialize variables
[n,D]       = size(Xnew) ;
Opt.GP      = OptGP ;
if Opt.GP.cov ~= inf
    d           = Opt.GP.cov;    % type of Martern
else
    d           = 1;
end
h1          = Opt.GP.h1 ;    % number of hyperparameters from covariance
h2          = Opt.GP.h2;     % number of hyperparameters from likelihood
hyp.cov     = hypVar(1:h1);
hyp.lik     = hypVar(h1+1:h1+h2);
ell         = exp(hyp.cov(1:D));
sf2         = exp(2*hypVar(D+1));
K           = zeros(n,n) ;

%% Calculate covariance matrix
for i = 1:D
    K = K_M(:,(i-1)*n+1:i*n) * d/ell(i)^2 + K;
end

if Opt.GP.cov ~= inf
    sqrtK = sqrt(K) ;
    expnK = exp(-sqrtK) ;
else
    expnK = exp(-1/2*K);
    sqrtK = [];
end

if      Opt.GP.cov == 3, t = sqrtK ; m =  (1 + t).*expnK;
elseif  Opt.GP.cov == 1,             m =  expnK;
elseif  Opt.GP.cov == 5, t = sqrtK ; m =  (1 + t.*(1+t/3)).*expnK;
elseif  Opt.GP.cov == inf,           m  = expnK;
end
K = sf2*m;
K =  K + eye(n)*exp(hyp.lik*2) ;
K = (K+K')/2 ; % This guarantees a symmetric matrix

%% Calculate inverse of covariance matrix
try
    CH = chol(K) ;
    invK = CH\(CH'\eye(n));
catch
    CH = chol(K+eye(n)*1e-4);
    invK = CH\(CH'\eye(n));
    warning('Covariance matrix in Nlikelihood is not positive semi-definite')
end

%% Calculate determinant of covariance matrix
logDetK = 2*sum(log(abs(diag(CH)))) ;

%% Calculate hyperperpriors
logprior = 0 ;
dlogpriorcov = zeros(1,h1) ;
for i = 1 : h1
    [A, dlogpriorcov(i)] =  priorGauss(Opt.GP.priorcov(1),Opt.GP.priorcov(2), hyp.cov(i) ) ;
    logprior = logprior + A ;
end

dlogpriorlik = zeros(1,h2);
for i = 1 : h2
    [A, dlogpriorlik(i)] =  priorGauss(Opt.GP.priorlik(1),Opt.GP.priorlik(2), hyp.lik(i) ) ;
    logprior = logprior + A ;
end

%% Calculate negative log-likeliehood
NLL = n/2*log(2*pi) + 1/2 * logDetK + 1/2 * Ynew'*invK*Ynew - logprior ;

%% Gradient calculation
if nargout == 2 % do only if No of output variable is 2 (if necessary)
    dsq_M = zeros(n,n*D) ;
    for i = 1 : D
        dsq_M(:,(i-1)*n+1:i*n)  = K_M(:,(i-1)*n+1:i*n) * (d)/ell(i)^2 ;
    end
    
    c = invK*Ynew;
    for i = 1 : h1
        dK = covMaternanisotropic(Opt.GP.cov,hyp.cov, sqrtK, expnK, dsq_M, Xnew, [], i);
        b = invK* dK ;
        dNLL_f.cov(i) = 1/2*trace(b) - 1/2*Ynew'*b*c ;
    end
    
    for i = 1 : h2
        dK = 2 * exp(hyp.lik(i)) * eye(n) * exp(hyp.lik(i));
        b = invK* dK ;
        dNLL_f.lik(i) = 1/2*trace(b) - 1/2*Ynew'*b*c ;
    end
    
    dNLL = [dNLL_f.cov';dNLL_f.lik'] - [dlogpriorcov';dlogpriorlik'];
end

return

function f = posterior_sample(Xnew,Ynew,Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2018-08-07.

% extration of variables from problem structure
nSpectralpoints = Opt.nSpectralpoints;
[n,D] = size(Xnew);
ell = exp(Opt.hyp.cov(1:D));
sf2 = exp(2*Opt.hyp.cov(D+1));
sn2 = exp(2*Opt.hyp.lik);

% Sampling of W and b
sW1  = lhsdesign(nSpectralpoints,D,'criterion','none');
sW2  = lhsdesign(nSpectralpoints,D,'criterion','none');
if Opt.cov ~= inf
    W = repmat(1./(ell)', nSpectralpoints, 1).*norminv(sW1).*sqrt(Opt.cov./chi2inv(sW2,Opt.cov));
else
    W = randn(nSpectralpoints,D) .* repmat(1./ell', nSpectralpoints, 1);
end

b = 2*pi*lhsdesign(nSpectralpoints,1,'criterion','none');

% Calculation of phi
phi = sqrt(2 * sf2 / nSpectralpoints) * cos(W * Xnew' + repmat(b, 1, n));

% Sampling of theta according to phi
A = phi * phi' + sn2 * eye(nSpectralpoints);

invA      = invChol(A);
mu_theta  = invA*phi*Ynew;
cov_theta = sn2*invA;
cov_theta = (cov_theta+cov_theta')/2;
theta     = mvnrnd(mu_theta,cov_theta)';

% Posterior sample (function) according to theta
f = @(x) (theta' * sqrt(2 * sf2 / nSpectralpoints) * cos(W * x' + repmat(b,1,size(x,1))))';
return

function f = mean_sample(Xnew,Ynew,Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2018-08-07.

% extration of variables from problem structure
nSpectralpoints = Opt.nSpectralpoints;
[n,D] = size(Xnew);
ell = exp(Opt.hyp.cov(1:D));
sf2 = exp(2*Opt.hyp.cov(D+1));
sn2 = exp(2*Opt.hyp.lik);

% Sampling of W and b
sW1  = lhsdesign(nSpectralpoints,D,'criterion','none');
sW2  = lhsdesign(nSpectralpoints,D,'criterion','none');
if Opt.cov ~= inf
    W = repmat(1./(ell)', nSpectralpoints, 1).*norminv(sW1).*sqrt(Opt.cov./chi2inv(sW2,Opt.cov));
else
    W = randn(nSpectralpoints,D) .* repmat(1./ell', nSpectralpoints, 1);
end

b = 2*pi*lhsdesign(nSpectralpoints,1,'criterion','none');

% Calculation of phi
phi = sqrt(2 * sf2 / nSpectralpoints) * cos(W * Xnew' + repmat(b, 1, n));

% Sampling of theta according to phi
A = phi * phi' + sn2 * eye(nSpectralpoints);

invA      = invChol(A);
mu_theta  = invA*phi*Ynew;

% Mean approximation (function) according to theta
f = @(x) (mu_theta' * sqrt(2 * sf2 / nSpectralpoints) * cos(W * x' + repmat(b,1,size(x,1))))';
return

function v = hypervolumemonte(P,r,N)
% Copyright (c) 2009, Yi Cao
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% HYPERVOUME    Hypervolume indicator as a measure of Pareto front estimate.
%   V = HYPERVOLUME(P,R,N) returns an estimation of the hypervoulme (in
%   percentage) dominated by the approximated Pareto front set P (n by d)
%   and bounded by the reference point R (1 by d). The estimation is doen
%   through N (default is 1000) uniformly distributed random points within
%   the bounded hyper-cuboid.
%
%   V = HYPERVOLUMN(P,R,C) uses the test points specified in C (N by d).
%
% See also: paretofront, paretoGroup

% Version 1.0 by Yi Cao at Cranfield University on 20 April 2008

% Example
%{
% an random exmaple
F=(randn(100,3)+5).^2;
% upper bound of the data set
r=max(F);
% Approximation of Pareto set
P=paretofront(F);
% Hypervolume
v=hypervolume(F(P,:),r,100000);
%}
% https://se.mathworks.com/matlabcentral/fileexchange/19651-hypervolume-indicator

% Check input and output
error(nargchk(2,3,nargin));
error(nargoutchk(0,1,nargout));

P=P*diag(1./r);
[n,d]=size(P);
if nargin<3
    N=1000;
end
if ~isscalar(N)
    C=N;
    N=size(C,1);
else
    C=rand(N,d);
end

fDominated=false(N,1);
lB=min(P);
fcheck=all(bsxfun(@gt, C, lB),2);

for k=1:n
    if any(fcheck)
        f=all(bsxfun(@gt, C(fcheck,:), P(k,:)),2);
        fDominated(fcheck)=f;
        fcheck(fcheck)=~f;
    end
end

v=sum(fDominated)/N;
return

function K = covMaternanisotropic(d, hyp, sqrtK,expnK, dsq_M, x, z, i)
% Copyright (c) 2005-2017 Carl Edward Rasmussen & Hannes Nickisch. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%    1. Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%    2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY CARL EDWARD RASMUSSEN & HANNES NICKISCH ``AS IS''
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL CARL EDWARD RASMUSSEN & HANNES NICKISCH OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The views and conclusions contained in the software and documentation
% are those of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of Carl Edward Rasmussen & Hannes Nickisch.
%
% The code and associated documentation is available from http://gaussianprocess.org/gpml/code.</pre>

[n,D] = size(x);
sf2 = exp(2*hyp(D+1));

if nargin<7                                                        % covariances
    if      d == 3, t = sqrtK ; m =  (1 + t).*expnK;
    elseif  d == 1,             m =  expnK;
    elseif  d == 5, t = sqrtK ; m =  (1 + t.*(1+t/3)).*expnK;
    elseif  d == inf, m = expnK;
    end
    K = sf2*m;
else                                                               % derivatives
    if i<=D                                               % length scale parameter
        Ki = dsq_M(:,(i-1)*n+1:i*n) ;
        if     d == 3,             dm = expnK;
        elseif d == 1, t = sqrtK ; dm = (1./t).*expnK;
        elseif d == 5, t = sqrtK ; dm = ((1+t)/3).*expnK;
        elseif d == inf; dm = -1/2*expnK;
        end
        
        K = sf2*dm.*Ki;
        K(Ki<1e-12) = 0;                                    % fix limit case for d=1
    elseif i==D+1                                            % magnitude parameter
        if      d == 3, t = sqrtK ; m =  (1 + t).*expnK;
        elseif  d == 1,             m =  expnK;
        elseif  d == 5, t = sqrtK ; m =  (1 + t.*(1+t/3)).*expnK;
        elseif  d == inf,           m = expnK;
        end
        K = 2*sf2*m;
    end
end
return

function D = sqdist(X1, X2)
% Copyright (c) 2016, Mo Chen
% All rights reserved.
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% * Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% * Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in
% the documentation and/or other materials provided with the distribution
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

% https://se.mathworks.com/matlabcentral/fileexchange/24599-pairwise-distance-matrix
% Pairwise square Euclidean distance between two sample sets
% Input:
%   X1, X2: dxn1 dxn2 sample matrices
% Output:
%   D: n1 x n2 square Euclidean distance matrix
% Written by Mo Chen (sth4nth@gmail.com).

D = bsxfun(@plus,dot(X2,X2,1),dot(X1,X1,1)')-2*(X1'*X2);
D(D<0) = 0 ; % check due to numerical errors
return

function [lp,dlp] = priorGauss(mu,s2,x)
% Copyright (c) 2005-2017 Carl Edward Rasmussen & Hannes Nickisch. All rights reserved.
%
% Redistribution and use in source and binary forms, with or without modification,
% are permitted provided that the following conditions are met:
%    1. Redistributions of source code must retain the above copyright notice,
%       this list of conditions and the following disclaimer.
%    2. Redistributions in binary form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY CARL EDWARD RASMUSSEN & HANNES NICKISCH ``AS IS''
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
% IN NO EVENT SHALL CARL EDWARD RASMUSSEN & HANNES NICKISCH OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
% EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%
% The views and conclusions contained in the software and documentation
% are those of the authors and should not be interpreted as representing official policies,
% either expressed or implied, of Carl Edward Rasmussen & Hannes Nickisch.
%
% The code and associated documentation is available from http://gaussianprocess.org/gpml/code.</pre>

lp  = -(x-mu).^2/(2*s2) - log(2*pi*s2)/2;
dlp = -(x-mu)/s2;
return

function [f,const] = pareto_objective(x,OptSample)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

Opt.Sample = OptSample;
f = zeros(size(x,1),size(Opt.Sample,2));
for i = 1:size(Opt.Sample,2)
    f(:,i) = Opt.Sample(i).f(x);
end
const = [];
return

function [Sample_pareto,Sample_xpareto,Sample_nadir] = Find_sample_pareto(Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

D = Opt.Gen.NoOfInputDim;               % Number of input dimensions
options = nsgaopt();                    % create default options structure
options.popsize = Opt.pop;              % populaion size
options.maxGen  = Opt.Generation;       % max generation
options.numObj = Opt.Gen.NoOfGPs;       % number of objectives
options.numVar = D;                     % number of design variables
options.numCons = 0;                    % number of constraints
options.outputfuns = [];                % saving pop
options.lb = zeros(1,D);                % lower bound of x
options.ub = ones(1,D);                 % upper bound of x
options.objfun = @(x) pareto_objective(x,Opt.Sample);        % objective function handle
options.useParallel = 'no';             % parallel computation is non-essential here
[~,result] = evalc('nsga2(options);');  % begin the optimization!

Sample_xpareto = zeros(Opt.pop,D);
Sample_pareto = zeros(Opt.pop,Opt.Gen.NoOfGPs);
result = result.pops(Opt.Generation,:);

for k = 1:Opt.pop
    Sample_xpareto(k,:) = result(k).var;
    Sample_pareto(k,:) = result(k).obj;
end

for k = 1:Opt.Gen.NoOfGPs
    Sample_nadir(k) = max(Sample_pareto(:,k));
end
return

function [Mean_pareto,Mean_xpareto] = Find_mean_pareto(Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

D = Opt.Gen.NoOfInputDim;               % Number of input dimensions
options = nsgaopt();                    % create default options structure
options.popsize = Opt.pop;              % populaion size
options.maxGen  = Opt.Generation;       % max generation
options.numObj = Opt.Gen.NoOfGPs;       % number of objectives
options.numVar = D;                     % number of design variables
options.numCons = 0;                    % number of constraints
options.outputfuns = [];                % saving pop
options.lb = zeros(1,D);                % lower bound of x
options.ub = ones(1,D);                 % upper bound of x
options.objfun = @(x) pareto_objective(x,Opt.Mean); % objective function handle
options.useParallel = 'no';             % parallel computation is non-essential here
[~,result] = evalc('nsga2(options);');  % begin the optimization!

Mean_xpareto = zeros(Opt.pop,D);
Mean_pareto  = zeros(Opt.pop,Opt.Gen.NoOfGPs);
result       = result.pops(Opt.Generation,:);

for k = 1:Opt.pop
    Mean_xpareto(k,:) = result(k).var;
    Mean_pareto(k,:) = result(k).obj;
end

return

function hv = hypervolume_2D(Yfront,r)
% The mex-file used follows the description of:
% 'M. Emmerich, K. Yang, A. Deutz, H. Wang and C. M. Fonseca. A
% Multicriteria Generalization of Bayesian Global Optimization.'
% Group website: http://liacs.leidenuniv.nl/~csmoda/index.php?page=code

AYfront = remove_points_above_reference(Yfront,r);
if isempty(AYfront)
    hv = 0;
else
    normvec = min(AYfront,[],1);
    A =  AYfront-repmat(normvec,size(AYfront,1),1);
    A =  A * diag(1./(r-normvec));
    A = -A + ones(size(A));
    A = sortrows(A,2);
    hyp_percentage = hypervolume2D(A,[0,0]);
    hv = prod(r-normvec)*hyp_percentage;
end
return

function hv = hypervolume_3D(Yfront,r)
% The mex-file used follows the description of:
% 'K. Yang, M. Emmerich, A. Deutz and C. M. Fonseca. A
% Computing 3-D Expected Hypervolume Improvement and Related Integrals in
% Asymptotically Optimal Time.'
% Group website: http://liacs.leidenuniv.nl/~csmoda/index.php?page=code

AYfront = remove_points_above_reference(Yfront,r);
if isempty(AYfront)
    hv = 0;
else
    normvec = min(AYfront,[],1);
    A = AYfront-repmat(normvec,size(AYfront,1),1) ;
    A = A *diag(1./(r-normvec));
    A = -A + ones(size(A)) ;
    A = sortrows(A,3);
    hyp_percentage = hypervolume3D(A,[0,0,0],[1,1,1]);
    hv = prod(r-normvec)*hyp_percentage;
end
return

function [index,hv_imp] = hypervolume_improvement_index(Ynew,Sample_nadir,Sample_pareto,Opt)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

r = Sample_nadir + 0.01*(max(Sample_pareto)-min(Sample_pareto));
index = [];

for i = 1 : Opt.NoOfBachSequential
    
    Yfront = Ynew(paretofront(Ynew),:);
    
    if size(Ynew,2) == 2
        hvY = hypervolume_2D(Yfront,r);
        
        for k = 1:size(Sample_pareto,1)
            A = [Ynew;Sample_pareto(k,:)];
            Afront = A(paretofront(A),:);
            hv = hypervolume_2D(Afront,r);
            hv_improvement(k) = hv-hvY;
        end
        
    elseif size(Ynew,2) == 3
        hvY = hypervolume_3D(Yfront,r);
        
        for k = 1:size(Sample_pareto,1)
            A = [Ynew;Sample_pareto(k,:)];
            Afront = A(paretofront(A),:);
            hv = hypervolume_3D(Afront,r);
            hv_improvement(k) = hv-hvY;
        end
        
    else
        AYfront = remove_points_above_reference(Yfront,r);
        normvec = min(AYfront,[],1);
        hyp_percentage = hypervolumemonte(AYfront-repmat(normvec,size(AYfront,1),1),r-normvec,3000);
        hvY = prod(r-normvec)*hyp_percentage;
        
        hv_improvement = zeros(size(Sample_pareto,1),1);
        for k = 1:size(Sample_pareto,1)
            B = [Ynew;Sample_pareto(k,:)];
            Bfront = B(paretofront(B),:);
            ABfront = remove_points_above_reference(Bfront,r);
            if isempty(ABfront)
                hv_improvement(k) = 0;
            else
                normvec = min(ABfront,[],1);
                hyp_percentage = hypervolumemonte(ABfront-repmat(normvec,size(ABfront,1),1),r-normvec,10000);
                hv = prod(r-normvec)*hyp_percentage;
                hv_improvement(k) = hv-hvY;
            end
        end
    end
    
    if i == 1
        hvY0 = hvY;
    end
    
    [~,Currentindex] = max(hv_improvement);
    Ynew = [Ynew;Sample_pareto(Currentindex,:)];
    index = [index;Currentindex];
end
hv_imp = hv_improvement(index(end))+hvY-hvY0;
return

function A = remove_points_above_reference(Afront,r)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

[A,~] = sortrows(Afront);
for p = 1:size(Afront,2);
    A = A(A(:,p)<=r(p),:);
end
return

function Opt = set_option_structure(Opt,X,Y)
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2017-13-12.

%% Extraction of input and output dimensions from X and Y
Opt.Gen.NoOfGPs = size(Y,2);           % Number of GPs to be trained
Opt.Gen.NoOfInputDim = size(X,2);      % Number of inputs dimensions
for i = 1 : Opt.Gen.NoOfGPs
    %% Set up GP options
    Opt.GP(i).cov = Opt.GP(i).matern;      % Matern type 1 / 3 / 5 / inf
    
    %% Set hyperpriors (MAP)
    Opt.GP(i).noiselimit = 0;                      % Upper bound on noise
    Opt.GP(i).var        = 10;                     % Upper bound on signal variance
    Opt.GP(i).h1         = Opt.Gen.NoOfInputDim+1; % Number of hyperparameters from covariance
    Opt.GP(i).h2         = 1;                      % Number of hyperparameters from likelihood
    
    %% priorGauss (mean, var)
    Opt.GP(i).priorlik  = [-6 ,Opt.GP(i).var];
    Opt.GP(i).priorcov  = [ 0 ,Opt.GP(i).var];
    
    %% Initial values for hyperparameters
    Opt.GP(i).hyp.cov       = zeros(1, Opt.GP(i).h1);
    Opt.GP(i).hyp.lik       = log(1e-2);
end
return
