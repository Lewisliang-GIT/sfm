function [best_P,inliers_idx] = estimate_T(xn,X,R,err_threshold)


alpha = 0.99;
epsilon = 0.0001;
update_T = @(epsilon) ceil(log10(1-alpha)/(log10(1-epsilon^2)));
T = update_T(epsilon);

iteration_count = 0;
while iteration_count <= T && iteration_count < 10000
    iteration_count = iteration_count +1;
    random_index = randperm(size(xn,2),2);
    
    t=estimate_T_DLT(R,xn(:,random_index),X(:,random_index));
    P = [R t];
    detal = pflat(P * [X;ones(1,size(X,2))]) - xn;
    detal = detal.^2;
    errors = detal(1,:) + detal(2,:);

    inliners = errors< err_threshold^2;

    temp_inliners_idx = find(inliners);
    
    inliners_number = size(temp_inliners_idx,2);
    
    n_epsilon = inliners_number / size(xn,2);
    if n_epsilon > epsilon
        epsilon = n_epsilon;
        best_P = P;
        inliers_idx = temp_inliners_idx;
        T = update_T(epsilon);
        iteration_count = 0;
    end 
end  
end
