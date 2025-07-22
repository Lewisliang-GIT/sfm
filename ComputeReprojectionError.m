function [err, res] = ComputeReprojectionError(P1, Xj, x1j)

    x1_proj = P1 * Xj;


    x1_proj = x1_proj / x1_proj(3);
  


    res1 = x1_proj - x1j;

    res = [res1(1:2)];
    err = norm(res)^2;
end
