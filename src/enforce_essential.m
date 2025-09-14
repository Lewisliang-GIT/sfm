function E = enforce_essential(E)

    [U,S,V] = svd(E);
    if det(U*V') < 0
        V=-V;
    end 
    E= U*diag([1 1 0])*V';
end