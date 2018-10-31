%Startup script for TS-EMO Algorithm

disp ('executing TS-EMO startup script...')
mydir = fileparts (mfilename ('fullpath')); % where am I located
addpath (mydir)
folders = {'Test_functions', 'NGPM_v1.4', 'Mex_files/hypervolume', 'Mex_files/invchol', 'Mex_files/pareto front', 'Direct'}; 
for d = folders, addpath (fullfile (mydir, d{1})), end
