function hypervolume = hypervolume2D(F,ub)
% Copyright (c) 2011, Johannes
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
% 
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
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
%
% Method for ND objective function values as described in:
%
% 'M. Fleischer. The measure of Pareto Optima Applications to 
%  Multi-objective Metaheuristics. EMO 2003, LNCSS 2632
%  519-533, 2003.'
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011
%
% Efficient method for 2D objective function values
    F  = -F' + ones(size(F'));
	L  = sortrows(F',1)';
	l  = length(L(1,:)); ub = ub + ones(1,size(L,1));
	hypervolume = 0;
	for i = 1:l
		hypervolume = hypervolume + ((L(1,i) - ub(1)) * (L(2,i) - ub(2)));
        ub(2)       = L(2,i);
    end
end
