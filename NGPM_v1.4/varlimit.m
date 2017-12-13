function var = varlimit(var, lb, ub)
% Function: var = varlimit(var, lb, ub)
% Description: Limit the variables in [lb, ub].
%
%         LSSSSWC, NWPU
%    Revision: 1.0  Data: 2011-04-20
%*************************************************************************

numVar = length(var);
for i = 1:numVar
    if( var(i) < lb(i) )
        var(i) = lb(i);
    elseif( var(i) > ub(i) )
        var(i) = ub(i);
    end
end

