# Pareto Set

version 1.1 by Yi Cao

Find the pareto set from n points with k objectives

It is motivated by Gianluca Dorini's isParetoSetMember program. The new m-file version is much faster than the 
C version because of the more elegant algorithm. The efficiency is significantly improved in version 3. 
By implementing a new sorting scheme and recoding to reduce overhead, the code is even faster than the mex version, 
paretomember, where sorting is not adopted. However, it puzzles me that the performance of mex code does not affected by sorting. 
Follow the links bellow to download these two codes for comparison. 

For the original unmodified file go to https://se.mathworks.com/matlabcentral/fileexchange/15181-pareto-set. 
