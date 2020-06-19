%% Example on how to use TSEMO_V4
% Copyright (c) by Eric Bradford, Artur M. Schweidtmann and Alexei Lapkin, 2020-23-05.

%% Step 1: Set path to folder
% Matlab Home, Set Path, Add folder with subfolders, select "TS-EMO-master"
% folder

%% Step 2: Create function file with multiobjective function to be optimized
% In this example we aim to find the Pareto front of "vlmop2", e.g. see
% vlmop2 for format
f = @vlmop2; % example function in the folder "Test_functions"

%% Step 3: Specify problem
no_outputs = 2;               % number of objectives
no_inputs  = 2;               % number of decision variables
lb = -2*ones(1,2);            % define lower bound on decision variables, [lb1,lb2,...]
ub =  2*ones(1,2);            % define upper bound on decision variables, [ub1,ub2,...]
      
%% Step 4: Generate initial dataset
dataset_size = 5*no_inputs;             % initial dataset size
X = lhsdesign(dataset_size,no_inputs);  % Latin hypercube design
Y = zeros(dataset_size,no_outputs);     % corresponding matrix of response data
for k = 1:size(X,1)
    X(k,:) = X(k,:).*(ub-lb)+lb;        % adjustment of bounds
    Y(k,:) = f(X(k,:));                 % calculation of response data
end

opt = TSEMO_options;             % call options for solver, see TSEMO_options file to adjust
opt.maxeval = 40;                % number of function evaluations before termination
opt.NoOfBachSequential = 1;      % number of function evaluations per iteration
% Total number of iterations = opt.maxeval/opt.NoOfBachSequential

%% Step 5: Start algorithm to find Pareto front
[Xpareto,Ypareto,X,Y,XparetoGP,YparetoGP,YparetoGPstd,hypf] = TSEMO_V4(f,X,Y,lb,ub,opt);

% INPUTS
% f denotes the function to be optimized
% X and Y   are the initial datasets to create a surrogate model
% lb and ub are the lower and upper bound of the decision variables
% opt is the option structure of the algorithm

% OUTPUTS
%   Xpareto and Ypareto correspond to the current best Pareto set and Pareto
%   front respectively. 
%   X and Y are the complete dataset of the decision variables and 
%   the objectives respectively. 
%   XparetoGP and YparetoGP represent the Pareto set and Pareto front of the 
%   final Gaussian process model within the algorithm. It is recommended to
%   use these as final result for problems with measurement noise. 
%   YparetoGPstd denotes the standard deviations of the predictions of
%   the GP pareto front YparetoGP  
%   hypf represents the final hyperparameters found for analysis

% For each iteration the current iteration number is displayed, the
% predicted hypervolume improvement and the time taken.

% TS-EMO creates a log file named "TSEMO_log.txt" that contains all relevant information 
% over the entire algorithm run. 

%% Step 6: Visualise results
figure
hold on
plot(Y(1:dataset_size,1),Y(1:dataset_size,2),'.','MarkerSize',14)
plot(Y(dataset_size+1:end,1),Y(dataset_size+1:end,2),'x','MarkerSize',8,'LineWidth',2)
plot(Ypareto(:,1),Ypareto(:,2),'O','MarkerSize',8,'LineWidth',2)
% plot(YparetoGP(:,1),YparetoGP(:,2),'x','MarkerSize',8,'LineWidth',2)
errorbar(YparetoGP(:,1),YparetoGP(:,2),YparetoGPstd(:,2),YparetoGPstd(:,2),YparetoGPstd(:,1),YparetoGPstd(:,1),'O','MarkerSize',8,'LineWidth',0.5)
legend('Initial LHC','Algorithm','Pareto front','GP Pareto front','Location','Northeast')
title('Results TS-EMO algorithm')
xlabel('f_1')
ylabel('f_2')
