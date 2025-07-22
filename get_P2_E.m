function [P_2,X,P_2s, Xs, best_point_count] = get_P2_E(E,x_1_n,x_2_n)
P_1 = [diag([1 1 1]) [0 0 0]'];

best_P2_idx = 1;
best_point_count = 0;
Xs = cell(1,4);
P_2s = extract_P_from_E(E);
for Pi=1:size(P_2s,2)
    X = triangulate_3D_point_DLT(x_1_n,x_2_n,P_1,P_2s{Pi});
    X=pflat(X);
    Xs{Pi} = X;
    in_front_P1 = X(3,:)>0;
    X2 = P_2s{Pi} *X;
    in_front_P2=X2(3, :) > 0;
    points_in_front = sum(in_front_P1 & in_front_P2);    
    if  points_in_front > best_point_count
        best_point_count = points_in_front;
        best_P2_idx = Pi;
    end 
end

P_2 = P_2s{best_P2_idx};
X=Xs{best_P2_idx};
end
