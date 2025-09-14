function J = ComputeJacobian(P, X,xn)
    
  J=[];
    for i = 1:size(xn,2)
        J=[J;
            (-P(3,:)*X(:,i)*[1 0 0]+P(1,:)*X(:,i)*[0 0 1])/(P(3,:)*X(:,i)).^2; ...
            (-P(3,:)*X(:,i)*[0 1 0]+P(2,:)*X(:,i)*[0 0 1])/(P(3,:)*X(:,i)).^2; ...
           ];
    end
end
