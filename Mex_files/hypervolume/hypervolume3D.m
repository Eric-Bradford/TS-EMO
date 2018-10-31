function hypervolume = hypervolume3D(F,ub,lb)
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

	hypervolume = 0;
    F  = -F' + ones(size(F'));
	[M, l] = size(F); ub = ub + ones(1,M); ub = ub';

	% Remove the duplicates from F and compute Lebesque measure of L
	L = unique(F','rows')';

	while l >= 1
		if (l > 1)
			b = zeros(M,1);
			spawn_vector = repmat(L(:,1), 1, M);
			for i = 1:M
				% Bound b(i) is either the least upper bound of the i-th value of the
				% other points or it is the value of the absolute upper bound ub(i)
				difL = (L(i,2:end) - L(i,1));
				lub = find(difL > 0);
				if (length(lub) > 0)
					b(i) = min(L(i,lub+1));
				else
					b(i) = ub(i);
				end

				b(i) = min((difL > 0) .* L(i,2:end) + (difL <= 0) * ub(i));
				% Update i-th spawn vector
				spawn_vector(i,i) = b(i);
			end

			% Compute lop-Off volume and update lebesgue measure
			lov = prod(b - L(:,1));
			hypervolume = hypervolume + lov;

			% Remove L(:,1) from L
			L = L(:,2:end);

			% Add the spawn_vector to L, but first filter dominated
			% solutions and the solutions that touch the upper bounds
			% from the spawn_vector
			L = nd_filter(L, spawn_vector, ub);

		else
			lov = prod(ub - L(:,1));
			hypervolume = hypervolume + lov;
			L = [];
		end
		% Update l
		[M, l] = size(L);
    end
end

function L = nd_filter(L, spawn_vector, ub)
% Implementation of the filter routine as described in:
%
% 'M. Fleischer. The measure of Pareto Optima Applications 
%  to Multi-objective Metaheuristics. EMO 2003, LNCSS 2632
%  519-533, 2003.'
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	[M, l_L] = size(L);
	[M, l_sp] = size(spawn_vector);
	do_assign = zeros(1, l_sp);

	for i = 1 : l_sp

		% Find if the spawnvector hits the upper bound
		at_ub = false;
		for j = 1:M
			if (spawn_vector(j,i) == ub(j))
				at_ub = true;
				break;
			end
		end

		% For this if statement, the following would be more elegant
		% (replacing the loop above), but less efficient:
		% if (all(spawn_vector(:,i) ~= ub))
		if (at_ub == false)
			do_assign(i) = 1;
			for j = 1 : l_L
				if (weakly_dominates(L(:,j), spawn_vector(:,i)))
					do_assign(i) = 0;
					break;
				end
			end
		end
	end
	L = [spawn_vector(:,find(do_assign == 1)), L];
end

function d = weakly_dominates(fA, fB)
% [d] = weakly_dominates(fA, fB)
%
% Compares two solutions A and B given their objective function
% values fA and fB. Returns whether A weakly dominates B.
%
% Input:
% - fA					- The objective function values of solution A
% - fB					- The objective function values of solution B
%
% Output:
% - d					- d is 1 if fA dominates fB, otherwise d is 0 
%
% Author: Johannes W. Kruisselbrink
% Last modified: March 17, 2011

	% Elegant, but not very efficient
	%d = (all(fA <= fB) && any(fA < fB));

	% Not so elegant, but more efficient
	d = true;
	for i = 1:length(fA)
		if (fA(i) > fB(i))
			d = false;
			return
		end
	end
end
