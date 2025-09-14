function errors = compute_epipolar_errors(F, x1s, x2s)
    l = F * x1s; % Compute epipolar lines in the second image
    l = l ./ sqrt ( repmat ( l (1 ,:).^2 + l (2 ,:).^2 ,[3 1])); % Normalize the lines
    errors = abs(sum(l .* x2s)); % Compute distances
end
