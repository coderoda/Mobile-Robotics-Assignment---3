function pts3D = algebraicTriangulation(x1, x2, P1, P2)
s1 = size(x1,2);
pts3D = zeros(4,s1);
for i = 1:s1
    MATRIX = zeros(4,4);
    MATRIX(1,:) = (x1(1,i).*P1(3,:)) - P1(1,:);
    MATRIX(2,:) = (x1(2,i).*P1(3,:)) - P1(2,:);
    MATRIX(3,:) = (x2(1,i).*P2(3,:)) - P2(1,:);
    MATRIX(4,:) = (x2(2,i).*P2(3,:)) - P2(2,:);
    [U,D,V] = svd(MATRIX);
    %     disp(size(V));
    pts3D(:,i) = V(:,4);
end
end