function [norm_F,T1_final,T2_final,inliers_a,inliers_b] = estimateFundamentalMatrixRansac(correspondingpts1, correspondingpts2)
%%% Homogeneous Representaion
% s1 = size(corresponding_pts1,1);
correspondingpts1(:,3) = ones(size(correspondingpts1,1),1);
correspondingpts2(:,3) = ones(size(correspondingpts1,1),1);

%%% Threshold error
error_thresh = 0.0010;
norm_F = zeros(3,3);

n = 8;
inliers = 0;

%%% RANSAC Loop
for i = 1:30000
    %%% Random 8 Points
    ind = randi(size(correspondingpts1,1), [n,1]);
    %%% normalizing those 8 points
    [norm2Dpts1 , T1] = normalize2DPoints(correspondingpts1(ind,:));
    [norm2Dpts2 , T2] = normalize2DPoints(correspondingpts2(ind,:));
    %%% Estimating F from the selected 8 points
    F_Estimate = T2' * Normalized_Estimate_F_Matrix(norm2Dpts1,norm2Dpts2) * T1;
    %%% Checking error from computed F
    myMat = correspondingpts1';
    myProd = F_Estimate * myMat;
%     err = sum((corresponding_pts2 .* (F_Estimate * corresponding_pts1')'),2);
    err = sum((correspondingpts2 .* myProd'),2);
%     current_inliers = size( find(abs(err) <= error_thresh) , 1);
    if (size( find(abs(err) <= error_thresh) , 1) > inliers)
        inliers = size( find(abs(err) <= error_thresh) , 1);
        norm_F = F_Estimate;
        T1_final = T1;
        T2_final = T2;
    end
end
disp(inliers);
myProd1 = norm_F * myMat;
err = sum((correspondingpts2 .* myProd1'),2);
[~,I]  = sort(abs(err),'ascend');
inliers_a = correspondingpts1(I(1:250),1:2);
inliers_b = correspondingpts2(I(1:250),1:2);
end