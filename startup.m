%Startup Script for TS-EMO Algorithm

%Add function folders to path
folders = {'Test_functions', 'NGPM_v1.4', 'Mex_files/hypervolume', 'Mex_files/invchol', 'Mex_files/pareto front' 'Direct'};
for i = 1:length(folders)
    new_path = strcat(pwd, '/', folders{i}, '/');
    addpath(new_path);
end