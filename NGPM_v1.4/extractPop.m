function nextpop = extractPop(opt, combinepop)
% Function: nextpop = extractPop(opt, combinepop)
% Description: Extract the best n individuals in 'combinepop'(population
%   size is 2n).
%
%         LSSSSWC, NWPU
%    Revision: 1.1  Data: 2011-07-12
%*************************************************************************

popsize = length(combinepop) / 2;
nextpop = combinepop(1:popsize);    %just for initializing

rankVector = vertcat(combinepop.rank);

n = 0;          % individuals number of next population
rank = 1;       % current rank number
idx = find(rankVector == rank);
numInd = length(idx);       % number of individuals in current front
while( n + numInd <= popsize )
    nextpop( n+1 : n+numInd ) = combinepop( idx );
    
    n = n + numInd;
    rank = rank + 1;
    
    idx = find(rankVector == rank);
    numInd = length(idx);
end

% If the number of individuals in the next front plus the number of individuals 
% in the current front is greater than the population size, then select the
% best individuals by corwding distance(NSGA-II) or preference distance(R-NSGA-II).
if( n < popsize )
    if(~isempty(opt.refPoints))
        prefDistance = vertcat(combinepop(idx).prefDistance);
        prefDistance = [prefDistance, idx];
        prefDistance = sortrows( prefDistance, 1);
        idxSelect  = prefDistance( 1:popsize-n, 2);       % Select the individuals with smallest preference distance
        nextpop(n+1 : popsize) = combinepop(idxSelect);
    else
        distance = vertcat(combinepop(idx).distance);
        distance = [distance, idx];
        distance = flipud( sortrows( distance, 1) );      % Sort the individuals in descending order of crowding distance in the front.
        idxSelect  = distance( 1:popsize-n, 2);           % Select the (popsize-n) individuals with largest crowding distance.
        nextpop(n+1 : popsize) = combinepop(idxSelect);
    end
end








