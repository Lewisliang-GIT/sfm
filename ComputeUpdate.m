function delta_Xj = ComputeUpdate(r, J, mu)

    delta_Xj = -inv((J'*J+mu*eye(size(J,2))))*J'*r;
end
