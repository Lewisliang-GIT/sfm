function X = triangulate_3D_point_DLT( points1, points2,camera1, camera2)
    X = [];
    for i= 1: size(points1,2)

       
        A=  [camera1(1,:)-points1(1,i)*camera1(3,:);
        camera1(2,:)-points1(2,i)*camera1(3,:);
        camera2(1,:)-points2(1,i)*camera2(3,:) ;
        camera2(2,:)-points2(2,i)*camera2(3,:) ];


        [~, ~, V] = svd(A);
        X_ = V(:, end);
        X = [X X_(1:4,:)];
    end
end