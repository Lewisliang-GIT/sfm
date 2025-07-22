function [H] = estimate_H_DLT(xs,ys)
num_points = size(xs, 2);
A = zeros(2 * num_points, 9);
for i = 1:num_points
    x = xs(1, i);
    y = xs(2, i);
    w = xs(3, i);

    x_prime = ys(1, i);
    y_prime = ys(2, i);
    w_prime = ys(3, i);

    A(2 * i - 1, :) = [-x, -y, -w,  0,  0,  0, x_prime * x, x_prime * y, x_prime * w];
    A(2 * i, :)     = [ 0,  0,  0, -x, -y, -w, y_prime * x, y_prime * y, y_prime * w];
end

[~, ~, V] = svd(A);

h = V(:, end);

H = reshape(h, [3, 3])';

H = H / H(3, 3);
end