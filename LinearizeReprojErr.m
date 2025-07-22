function [r, J] = LinearizeReprojErr(P1, Xj, x1j)
   
    x1_proj = P1 * Xj;

    x1_proj = x1_proj / x1_proj(3);

  
    r1 = x1_proj - x1j;
    r = reshape(r1(1:2,:),[],1);

  
    J1 = ComputeJacobian(P1, Xj,x1j);

    J = J1;
end

