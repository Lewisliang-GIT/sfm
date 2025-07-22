function [best_E, inliers]= estimate_E_robust(x1, x2, eps)
alpha = 0.9999;
epsilon = 0.1;
T = log10(1-alpha)/(log10(1-epsilon^8));
ri = randperm(size(x1, 2), 8);
rp1 = x1(:, ri);
rp2 = x2(:, ri);
Eapp1 = estimate_F_DLT(rp1, rp2);
best_E = enforce_essential(Eapp1);
while T>0
    sample_idx = randperm(size(x1, 2), 8);
    x1_sample = x1(:, sample_idx);
    x2_sample = x2(:, sample_idx);
    E_candidate = enforce_essential ( estimate_F_DLT ( x1_sample , x2_sample));
    errors = (compute_epipolar_errors(E_candidate',x2,x1).^2 ...
        + compute_epipolar_errors(E_candidate,x1,x2).^2)/2;
    inliners = errors < eps^2;
    nInliers = sum(inliners);
    n_epsilon = nInliers / size(x1,2);
    if n_epsilon > epsilon
        epsilon = n_epsilon;
        best_E = E_candidate;
        inliers = find(inliners);
        T = log10(1-alpha)/(log10(1-epsilon^8));
    end
    T=T-1;
end
end