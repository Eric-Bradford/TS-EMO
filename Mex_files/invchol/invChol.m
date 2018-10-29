function invA = invChol(A)
    CH = chol(A);
    invA = CH\(CH'\eye(n));
end
    
