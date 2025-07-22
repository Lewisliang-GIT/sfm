function [R,T,best_X, inliers_idx] = parallel_RANSAC( ...
    x1,x2,e_threshold,h_threshold)

alpha = 0.999;

epsilon_e = 0.001;
epsilon_h = 0.001;

update_T = @(epsilon,sample_point_num) ceil(log10(1-alpha)/(log10(1-epsilon^sample_point_num)));
build_skew = @(t) [0 -t(3) t(2); t(3) 0 -t(1); -t(2) t(1) 0];

T_e = update_T(epsilon_e,8);
T_h = update_T(epsilon_h,4);

iteration_e = 0;
iteration_h = 0;
iteration=0;

update_flag = 0;
while iteration_e <= T_e || iteration_h <= T_h || iteration <1000
    iteration=iteration+1;
    if iteration_e <= T_e
        iteration_e = iteration_e +1;
        random_index = randperm(size(x1,2),8);

        [E_n] = estimate_F_DLT( ...
            x1(:,random_index),x2(:,random_index));
        bE_n = enforce_essential(E_n);

        errors_e = (compute_epipolar_errors(bE_n',x2,x1).^2 ...
            + compute_epipolar_errors(bE_n,x1,x2).^2)/2;

        inliners_e = errors_e< e_threshold^2;

        temp_inliners_idx_e = find(inliners_e);

        inliners_number_e = size(temp_inliners_idx_e,2);

        n_epsilon_e = inliners_number_e / size(x1,2);
        if n_epsilon_e > epsilon_e
            epsilon_e = n_epsilon_e;
            best_E = bE_n;
            inliers_idx = temp_inliners_idx_e;
            T_e = update_T(epsilon_e,8);

            [P2,best_X,~,~,~]= get_P2_E( ...
                best_E,x1(:,inliers_idx),x2(:,inliers_idx));
            iteration_e=0;
            R= P2(:,1:3);
            T= P2(:,4);
            update_flag = 0;
        end
    end

    % estimate H
    if iteration_h <= T_h
        iteration_h = iteration_h +1;
        random_index = randperm(size(x1,2),4);
        H=estimate_H_DLT(x1(:,random_index),x2(:,random_index));
        detal = pflat(H*x1) - x2;
        detal = detal.^2;
        errors = detal(1,:) + detal(2,:);

        inliners_h = errors< h_threshold^2;
        temp_inliners_idx_h = find(inliners_h);

        inliners_number_h = size(temp_inliners_idx_h,2);

        n_epsilon_h = inliners_number_h / size(x1,2);
        if n_epsilon_h > epsilon_h
            % R T from H
            [R1,t1,R2,t2]=homography_to_RT(H,x1,x2);
            E1 = build_skew(t1) * R1;
            E1 = enforce_essential(E1);
            if rank(E1)==2
                [P2_1,~,~,~,count_1]= get_P2_E( ...
                    E1,x1(:,temp_inliners_idx_h),x2(:,temp_inliners_idx_h));
            end
            E2 = build_skew(t2) * R2;
            E2 = enforce_essential(E2);
            if rank(E2)==2
                [P2_2,~,~,~,count_2]= get_P2_E( ...
                    E2,x1(:,temp_inliners_idx_h),x2(:,temp_inliners_idx_h));
            end

            %select best e
            best_h_E = E1;
            best_h_P = P2_1;
            if count_2 > count_1
                best_h_E = E2;
                best_h_P = P2_2;
            end

            R= best_h_P(:,1:3);
            T= best_h_P(:,4);

            errors_e_h = (compute_epipolar_errors(best_h_E',x2,x1).^2 ...
                + compute_epipolar_errors(best_h_E,x1,x2).^2)/2;

            inliners_e_h = errors_e_h< e_threshold^2;

            temp_inliners_idx_e_h = find(inliners_e_h);

            inliners_number_e_h = size(temp_inliners_idx_e_h,2);

            n_epsilon_e_h = inliners_number_e_h / size(x1,2);

            %             E from H is better
            if n_epsilon_e_h > epsilon_e
                epsilon_e = n_epsilon_e_h;
                best_E = best_h_E;
                inliers_idx = temp_inliners_idx_e_h;
                T_e = update_T(epsilon_e,8);
                iteration_e=0;
                [P2,best_X,~,~,~]= get_P2_E( ...
                    best_E,x1(:,inliers_idx),x2(:,inliers_idx));
                R= P2(:,1:3);
                T= P2(:,4);
                update_flag = 1;
            end

            epsilon_h = n_epsilon_h;
            T_h = update_T(epsilon_h,4);
            iteration_h=0;
        end
    end

end
if update_flag == 1
    fprintf("from H better\n");
else
    fprintf("from E better\n");
end


end % function end