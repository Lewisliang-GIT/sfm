function F = enforce_fundamental(F_approx)
    [U, D, V] = svd(F_approx);
    D(end, end) = 0; 
    F = U * D * V';
end
