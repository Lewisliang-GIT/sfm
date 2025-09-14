function P2_options = extract_P_from_E(E)
    % Decompose E
    [U, ~, V] = svd(E);
    if det(U * V') < 0
        V = -V;
    end
    
    % Define W
    W = [0 -1 0; 1 0 0; 0 0 1];
    
    % Compute camera matrices
    P2_options = {
        [U * W * V', U(:, 3)], ...
        [U * W * V', -U(:, 3)], ...
        [U * W' * V', U(:, 3)], ...
        [U * W' * V', -U(:, 3)]
    };
end
