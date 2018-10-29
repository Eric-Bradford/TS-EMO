function invA = invChol(A)
    n    = size(A,2);
    CH   = chol(A);
    invA = CH\(CH'\eye(n));
end
    
