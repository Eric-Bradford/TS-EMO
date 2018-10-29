function invA = invChol(A)
    CH = chol(K);
    invK = CH\(CH'\eye(n));
end
    