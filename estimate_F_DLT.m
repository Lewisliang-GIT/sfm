function F = estimate_F_DLT(x1, x2)

    n_points = size(x1, 2);
    M = zeros(n_points, 9);
    for i = 1:n_points
        xx = x2(:, i) * x1(:, i)';
        M(i, :) = xx(:)';
    end

    [U, S, V] = svd(M);
    v = V(:, end); 
    F = reshape(v, [3, 3]);

    nMv =  norm(M*v);
    sigluars = diag(S)';

    
end
