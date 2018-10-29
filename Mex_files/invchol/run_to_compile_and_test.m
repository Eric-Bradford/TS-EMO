% Description:
%       Computes the inverse of a symmetric matrix using the cholesky
%       decomposition and compares the runtime to that of the standard inv
%       function.
%
% Author: 
%		Eric Blake
%		Black River Systems Company
%		blake@brsc.com
%		05/14/2012


% Compile the MEX file.
clc;
disp('Compiling invChol_mex.c...');

if strcmpi(computer, 'PCWIN64') || strcmpi(computer, 'GLNXA64') || ...
        strcmpi(computer, 'MACI64')
    disp('64-bit detected, compiling with -largeArrayDims flag...');
    mex ('invChol_mex.c', ...
        ['-I' matlabroot '/extern', '/examples', '/refbook'],...
        [matlabroot '/extern', '/examples', '/refbook', '/fort.c'], ...
        '-lmwlapack', '-largeArrayDims');
else
    disp('32-bit detected, compiling without -largeArrayDims flag...');
    mex('invChol_mex.c',...
        ['-I' matlabroot '/extern', '/examples', '/refbook'],...
        [matlabroot '/extern', '/examples', '/refbook', '/fort.c'], ...
        '-lmwlapack');
end

% Number of multiplies to test.
n_mults = 10;

% Size of matrix to test.
matdim = 1000;

% Generate random symmetric matrix.
if verLessThan('matlab', '7.12')
    rng(0);  % Set random seed.
else
    rand('seed', 0);  % Set random seed in pre R2011a version.
end
disp(['Generating random ' num2str(matdim) ' x ' num2str(matdim) ...
    ' symmetric matrix....']);
A = rand(matdim, matdim);
A = A'*A;

% Compute the inverse n_mults times using standard inv functin (uses LU).
disp('Executing traditional MATLAB inv...');
drawnow; % Ensures previous line is executed.
tic;
for i=1:n_mults
    inv(A);
end
t_standard=toc;
disp(['Runtime: ' num2str(t_standard) ' seconds']);

% Compute the inverse n_mults times using cholesky inverse.
A_cholinv = invChol_mex(A); % Run one time first to catch any bugs.
disp('Executing invChol_mex...');
drawnow; % Ensures previous line is executed.
tic;

for i=1:n_mults
    invChol_mex(A);
end
t_chol=toc;
disp(['Runtime: ' num2str(t_chol) ' seconds']);

% Compute and display the max numerical difference.
err = A_cholinv - inv(A);
str =sprintf('Maximum numerical difference between inv and invChol_mex: %d', ...
    max((err(:))));
disp(str);




