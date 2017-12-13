function opt = verifyOpt(opt)
% Function: opt = verifyOpt(opt)
% Description: Verify the optimization model.
%         LSSSSWC, NWPU
%   Revision: 1.1  Data: 2011-07-15
%*************************************************************************


%*************************************************************************
%  popsize
%*************************************************************************
if ( mod(opt.popsize, 2) ~= 0 )
    warning('NSGA2:PopSizeError', 'The population size shoud be even number!%d => %d', opt.popsize, opt.popsize+1);
    opt.popsize = opt.popsize + 1;
end

%*************************************************************************
% lb, ub
%*************************************************************************
if( length(opt.lb)~=opt.numVar || length(opt.lb)~=opt.numVar )
    error('NSGA2:OptModelError', 'The numbers of lower and upper bounds(%d,%d) should be equal to the design variable number(%d)!', ...
        length(opt.ub), length(opt.lb), opt.numVar);
end

%*************************************************************************
% vartype
%*************************************************************************
if( length(opt.vartype) ~= opt.numVar )
    warning('NSGA2:OptModelWarning', 'Design variables'' data type error! All the type is set to REAL coding (vartype=1)!');
    opt.vartype = ones(1, opt.numVar);
end

%*************************************************************************
% nameObj, nameVar, nameCons
%*************************************************************************
if( ~iscell(opt.nameObj) || ~iscell(opt.nameVar) || ~iscell(opt.nameCons))
    error('NSGA2:OptModelError', 'The names of objectives, design variables or constraints should be specified in cell array, for example, {''obj1'',''obj2''}');
end

if( (~isempty(opt.nameObj)  && length(opt.nameObj)~=opt.numObj) || ...
    (~isempty(opt.nameVar)  && length(opt.nameVar)~=opt.numVar) || ...
    (~isempty(opt.nameCons) && length(opt.nameCons)~=opt.numCons))
    error('NSGA2:OptModelError', 'All names of objectives, design variables or constraints should be specified, if one is specified!');
end

%*************************************************************************
% useparallel
%*************************************************************************
if( ~ischar(opt.useParallel) || ...
    isempty( find(strcmpi(opt.useParallel, {'yes', 'no'}))) )
    error('NSGA2:OptParamError', 'useParallel can be only "yes" or "no"!');
end

%*************************************************************************
% R-NSGA-II parameters
%*************************************************************************
% refPoints
if( ~isempty(opt.refPoints) && size(opt.refPoints,2)~=opt.numObj)
    error('NSGA2:OptParamError', 'The reference points has the format refPoints(nPoint, numObj)!');
end
% refWeight
if( ~isempty(opt.refPoints) && ~isempty(opt.refWeight) && length(opt.refWeight)~=opt.numObj)
    error('NSGA2:OptParamError', 'The weight factor vector used in R-NSGA-II must has the length of numObj!');
end









