
function run_sfm(id)
rng(42);
addpath Z:/Downloads/vlfeat-0.9.21/toolbox/;
vl_setup()

[K, img_names, init_pair, pixel_threshold] = get_dataset_info(id);
epipolar_threshold = pixel_threshold/K(1,1);
homography_threshold = 3 * pixel_threshold/K(1,1);
translation_threshold = 3 * pixel_threshold/K(1,1);
fs=cell(1,size(img_names,2));
ds=cell(1,size(img_names,2));
for i=1:size(img_names,2)
    image_i = imread(img_names{i});
    [fi, di] = vl_sift(single(rgb2gray(image_i)), 'PeakThresh', 1);
    fs{i} = fi;
    ds{i} = di;
end

%step 1
R_s = cell(1,size(img_names,2));
R_s{1} = [1 0 0; 0 1 0; 0 0 1];
for i = 1 : size(img_names,2)-1

    f1 = fs{i};
    f2 = fs{i+1};
    d1 = ds{i};
    d2 = ds{i+1};
    [matches, ~] = vl_ubcmatch(d1,d2);
    x1 = [f1(1,matches(1,:)); ...
        f1(2,matches(1,:)); ...
        ones(1,size(f1(2,matches(1,:)),2))];
    x2 = [f2(1,matches(2,:)); ...
        f2(2,matches(2,:)); ...
        ones(1,size(f2(2,matches(2,:)),2))];
    x1_normalized = inv(K)*x1;
    x2_normalized = inv(K)*x2;

    %stander 8-point algorithm
    % [E, inliers_idx] = estimate_E_robust(x1_normalized,x2_normalized,epipolar_threshold);
    % [P2,X,~,~]= get_P2_E(E,x1_normalized(:,inliers_idx),x2_normalized(:,inliers_idx));
    % R_s{i+1} = P2(:,1:3);

    %parallel RANSAC
    [R,~,~,~] = parallel_RANSAC( ...
        x1_normalized,x2_normalized, ...
        epipolar_threshold,homography_threshold);
    R_s{i+1} = R;

end

%step 2
R_s_abs = cell(1, size(R_s,2));
R_s_abs{1} = [1 0 0; 0 1 0; 0 0 1];
for i = 2 : size(R_s,2)
    R_s_abs{i} = R_s{i} * R_s_abs{i-1};
end

%step 3
%load initial pair
image_1 = imread(img_names{init_pair(1)});
image_2 = imread(img_names{init_pair(2)});
f1 = fs{init_pair(1)};
d1 = ds{init_pair(1)};
f2 = fs{init_pair(2)};
d2 = ds{init_pair(2)};

[matches, scores] = vl_ubcmatch(d1,d2);
x1 = [f1(1,matches(1,:)); f1(2,matches(1,:)); ones(1,size(f1(2,matches(1,:)),2))];
x2 = [f2(1,matches(2,:)); f2(2,matches(2,:)); ones(1,size(f2(2,matches(2,:)),2))];

x1_normalized = inv(K)*x1;
x2_normalized = inv(K)*x2;


%stander 8-point algorithm
% [E, inliers_idx] = estimate_E_robust(x1_normalized,x2_normalized,epipolar_threshold);
% [P2,X,~,~,~]= get_P2_E(E,x1_normalized(:,inliers_idx),x2_normalized(:,inliers_idx));


%parallel RANSAC
[R,T,X,inliers_idx] = parallel_RANSAC( ...
    x1_normalized,x2_normalized, ...
    epipolar_threshold,homography_threshold);

X_0 = R_s_abs{init_pair(1)}'*X(1:3,:);

d1_matches = d1(:,matches(1,:));
desc_x = x1(:,inliers_idx);
desc_X = d1_matches(:,inliers_idx);
desc_img = image_1;

%step 4
Ps = cell(1,size(img_names,2));
xs = cell(1,size(img_names,2));
Xs = cell(1,size(img_names,2));
for i=1:size(img_names,2)
    image_i = imread(img_names{i});
    fi = fs{i};
    di = ds{i};
    [matches, scores] = vl_ubcmatch(di,desc_X);
    xi = [fi(1,matches(1,:)); fi(2,matches(1,:)); ones(1,size(fi(2,matches(1,:)),2))];
    xi = inv(K)*xi;
    Xi =X_0(:,matches(2,:));

    [Ps{i}, inliers_t_idx] = estimate_T(xi,Xi,R_s_abs{i},translation_threshold);
    Xs{i} = Xi(:,inliers_t_idx);
    xs{i} = xi(:,inliers_t_idx);

end


%step 5
for i=1:size(Ps,2)
    [P] = refine_T(Ps{i},xs{i},[Xs{i};ones(1,size(Xs{i},2))]);
    Ps{i} = P;
end

%step 6
Xs= cell(1,size(img_names,2)-1);
for i=1:size(img_names,2)-1
    image_1 = imread(img_names{i});
    image_2 = imread(img_names{i+1});
    f1 = fs{i};
    d1 = ds{i};
    f2 = fs{i+1};
    d2 = ds{i+1};

    [matches, scores] = vl_ubcmatch(d1,d2);
    x1 = [f1(1,matches(1,:)); f1(2,matches(1,:)); ones(1,size(f1(2,matches(1,:)),2))];
    x2 = [f2(1,matches(2,:)); f2(2,matches(2,:)); ones(1,size(f2(2,matches(2,:)),2))];
    x1n = inv(K)*x1;
    x2n = inv(K)*x2;


    Xi=triangulate_3D_point_DLT(x1n,x2n,Ps{i},Ps{i+1});

    Xi=pflat(Xi);

    Xi_mean = mean(Xi')';
    Xi_dist = vecnorm(Xi - Xi_mean,2,1);
    threshold = 2 * quantile(Xi_dist,0.90);
    Xi_inliners = Xi_dist <= threshold;
    Xi = Xi(:,Xi_inliners);

    Xs{i} = Xi;
end

figure
for i=1: size(Xs,2)
    plot3(Xs{i}(1,:),Xs{i}(2,:),Xs{i}(3,:),'.', 'MarkerSize',3);
    hold on
end
plotcams(Ps);
hold off
axis equal

end





