function [Pb] = refine_T(P,xn,X)

    % R = P(:,1:3);
    % t = P(:,4);
    mu = 0.01;
    max_iters = 100000;
    
    
    num_points = size(X, 2);
    P_refined = P;
    err_before = zeros(1,size(X,2));
    for i = 1:num_points
        [err, ~] = ComputeReprojectionError(P_refined,  X(:, i), xn(:, i));
        err_before(i) = err;
    end
    t_err_before=sum(err_before);
    m_err_before = median(err_before);
    
        for iter = 1:max_iters
            [r, J] = LinearizeReprojErr(P_refined, X, xn);
          
            % Compute update
            delta_t =-pinv((J'*J+mu*eye(size(J,2))))*J'*r;
            [err, ~] = ComputeReprojectionError(P_refined, X,  xn);
            R = P_refined(:,1:3);
            
            t = P_refined(:,4)+delta_t;
            P_refined_delta=[R t];
            [err_re, ~] = ComputeReprojectionError(P_refined_delta, X, xn);
            if err_re< err
                P_refined =P_refined_delta;
                mu = mu/10;
            else
                mu = mu*10;
            end
            if norm(delta_t) < 1e-6
                break;
            end
        end
    
    err_after = 0;
    for i = 1:num_points
        [err, ~] = ComputeReprojectionError(P_refined,X(:, i), xn(:, i));
        err_after(i) = + err;
    end
    t_err_after=sum(err_after);
    m_err_after = median(err_after);
    Pb = P_refined;
    end