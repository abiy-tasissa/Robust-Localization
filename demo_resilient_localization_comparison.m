% ----------------------------------------------------------------------------------------------------
% This script is on the problem of localization, where the objective is 
% to determine the positions of target nodes given their distances to m 
% anchor nodes with known positions. We assume that the last anchor node 
% serves as the central node.
%
% In our setting, a random subset (alpha) of the target nodes is highly 
% corrupted. This means that, for these nodes, a random subset (k) of their distance 
% measurements to the anchors is highly corrupted. 
%
% We compare our algorithm with the SR-Hybrid algorithm in the ref. below:
% Zaeemzadeh, Alireza, et al. "Robust target localization based on squared range iterative reweighted least squares." 
% 2017 IEEE 14th International Conference on Mobile Ad Hoc and Sensor Systems (MASS). IEEE, 2017.
% ------------------------------------------------------------------------------------------------------
% Set number of anchors, embedding dimensions
% m = number of anchors
% r = embedding dimension
% alpha = total number of highly corrupted nodes
% k = number of corrupted measurements (per a corrupted node)
% beta1 = noise parameter for normal nodes
% beta2 = noise parameter for highly corrupted nodes
m = 15;
r = 2;
alpha = 4;
k = 3;
beta1_min = -100;
beta1_max= 100;
beta2_min = -100;
beta2_max= 100;
beta = 0.75;
% standadrd deviation of Gaussian noise
std_dev = 40;
% If central node is corrrupted or not
central_corrupt = 1;
% Number of near and far target nodes
opts.num_near = 50; 
opts.num_far = 50; 
% num_repeats (since there is randomness)
num_repeats= 50;
% Create arrays to store errors
mean_rel_err1 = zeros(num_repeats,1);
mean_rel_err2 = zeros(num_repeats,1);
mean_pos_err1 = zeros(num_repeats,1);
mean_pos_err2 = zeros(num_repeats,1);
mean_anchors_dist1 = zeros(num_repeats,1);
mean_anchors_dist2 = zeros(num_repeats,1);
mean_dist_err = zeros(num_repeats,1);
identification_error = zeros(num_repeats,1);
identification_error1 = zeros(num_repeats,1);
naive_error = zeros(num_repeats,1);
coherence = zeros(num_repeats,1);
for i = 1:num_repeats
    i
    % Generate the anchors, target nodes (both near and far)
    % near = target nodes that are near to the anchors
    % far = target nodes that are far to the anchors
    [anchors, target_near, target_far] = generate_near_far_nodes(r,m,opts);
    % Form squared distance matrix of target nodes
    nodes_all = [anchors target_near target_far];
    % Squared distance matrix
    dist = squareform(pdist(nodes_all'));
    D = dist.*dist;
    % Apply corruption on (1,2) block of the distance matrix
    sz_target_near = size(target_near);
    num_near_final = sz_target_near(2);
    sz_target_far = size(target_far);
    num_far_final = sz_target_far(2);
    % Too few far nodes: break the loop
    if num_far_final<alpha
        continue
    end
    F = D(1:m,m+1:end);
    % Pick k_2 indices at random to corrupt
    rand_idx = randperm(num_far_final);
    rand_idx = rand_idx(1:alpha);
    rand_idx_columns = zeros(k,alpha);
    % Corruption of normal nodes
    % Uniform noise
    beta_normal = beta1_min + (beta1_max - beta1_min) * rand(size(F));
    gau_noise= std_dev*randn(size(F));
    if central_corrupt==1
        sqrt_Fc = sqrt(F) + (1-beta)*gau_noise+beta*beta_normal;
        F_corrupted = sqrt_Fc.*sqrt_Fc;
    else
        sqrt_Fc = sqrt(F) + (1-beta)*gau_noise+beta*beta_normal;
        F_corrupted = sqrt_Fc.*sqrt_Fc;
        F_corrupted(m,:) = F(m,:);
    end   
    % Generate highly corrupted nodes
    for j = 1:alpha
        if central_corrupt ==1
            rand_idx_1 = randperm(m);
            rand_idx_1 = rand_idx_1(1:k);
        else
            rand_idx_1 = randperm(m-1);
            rand_idx_1 = rand_idx_1(1:k);
        end
        rand_idx_columns(:,j) = rand_idx_1';
        F1 = F(rand_idx_1,rand_idx(j)+num_near_final);
        beta_corrup = beta2_min + (beta2_max - beta2_min) * rand(size(F1));
        gau_noise = std_dev*randn(size(F1));
        sqrt_Fc = sqrt(F1) + (1-beta)*gau_noise+beta*beta_corrup; 
        F_corrupted(rand_idx_1,rand_idx(j)+num_near_final) = sqrt_Fc.*sqrt_Fc;
    end
    central = anchors(:,m);
    % Re-define anchors removing central node
    anchors = anchors(:,1:m-1);
    % Create the system Aq=b and Aq = b_corrupted
    A = (anchors-central)';
    anchor_norms = vecnorm(anchors).^2;
    % Apply sparse+dense algorithm to identified nodes, to remove
    % outliers
    rel_err1 = zeros(alpha,1);
    rel_err2 = zeros(alpha,1);
    anchors_dist1 = zeros(alpha,1);
    anchors_dist2 = zeros(alpha,1);
    pos_err1 = zeros(alpha,1);
    pos_err2 = zeros(alpha,1);
    dist_err = zeros(alpha,1);
    far_nodes = zeros(r,alpha);
    for tt = 1:alpha
        if central_corrupt==1
        % Distance from target to central node
            t_c_dist = F_corrupted(m,rand_idx(tt)+num_near_final);
        else
            t_c_dist = norm(target_far(:,rand_idx(tt))-central)^2;
        end
        b= 0.5*(t_c_dist - F_corrupted(1:m-1,rand_idx(tt)+num_near_final) + anchor_norms' -norm(central)^2);
        q = A\b;
        % Call algorithm 1 and algorithm 2
        R = sqrt(F_corrupted(1:m-1,rand_idx(tt)+num_near_final));
        [q_est1,~] = PUpositionSRHybrid(A(:,1),A(:,2),R,std_dev);
        q_est2= robust_alg(A,b) ;
        % Performance metrics
        target_pos= target_far(:,rand_idx(tt));
        pos_err1(tt) = norm(target_pos-q_est1)^2;
        pos_err2(tt) = norm(target_pos-q_est2)^2;
        rel_err1(tt) = norm(target_pos-q_est1)/norm(target_pos);
        rel_err2(tt) = norm(target_pos-q_est2)/norm(target_pos);
        t_a_vec1 = q_est1-[anchors central];
        t_a_vec2 = q_est2-[anchors central];
        tar_idx = rand_idx(tt)+num_near_final;
        dist_err(tt) = sum((sqrt(F(1:m,tar_idx))-sqrt(F_corrupted(1:m,tar_idx))).^2);
        anchors_dist1(tt) = mean(vecnorm(t_a_vec1))/mean(sqrt(F(1:m,tar_idx)));
        anchors_dist2(tt) = mean(vecnorm(t_a_vec2))/mean(sqrt(F(1:m,tar_idx)));         
    end
    mean_rel_err1(i) = mean(rel_err1); 
    mean_rel_err2(i) = mean(rel_err2);
    mean_pos_err1(i) = mean(pos_err1); 
    mean_pos_err2(i) = mean(pos_err2);
    mean_dist_err(i) = sum(dist_err)/(alpha*m);
    mean_anchors_dist1(i) = mean(anchors_dist1);
    mean_anchors_dist2(i) = mean(anchors_dist2);
    
end
