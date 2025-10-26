% --------------------------------------------------------------------------------------
% This script is on the problem of localization, where the objective is 
% to determine the positions of target nodes given their distances to m 
% anchor nodes with known positions. We assume that the last anchor node 
% serves as the central node.
%
% In our setting, a random subset (alpha) of the target nodes is highly 
% corrupted. This means that, for these nodes, a random subset (k) of their distance 
% measurements to the anchors is highly corrupted. 
%
% The corruption model follows:
%   (1 +  U(beta_min,beta_max)) * d_{i,j}^2
% where U(beta_min,beta_max) is uniform random variable in the interval [beta_min,beta_max], and 
% d_{i,j}^2 represents the squared distance. 
%
% Our goal is to:
%   (a) Identify the highly corrupted nodes.
%   (b) Estimate their positions in a robust manner.
% --------------------------------------------------------------------------------------
% Set number of anchors, embedding dimensions
% m = number of anchors
% r = embedding dimension
% alpha = total number of highly corrupted nodes
% k = number of corrupted measurements (per a corrupted node)
% beta1 = noise parameter for normal nodes
% beta2 = noise parameter for highly corrupted nodes
m = 9;
r = 2;
alpha = 4;
k = 3;
beta1_min = 0.00;
beta1_max= 0.05;
beta2_min= 0.15;
beta2_max = 0.20;
% If central node is corrrupted or not
central_corrupt = 1;
% Options to visualize
opts.visualize = 0;
% Number of near and far target nodes
opts.num_near = 50; 
opts.num_far = 50; 
% num_repeats (since there is randomness)
num_repeats= 50;
% Create arrays to store errors
mean_rel_err = zeros(num_repeats,1);
mean_pos_err = zeros(num_repeats,1);
mean_anchors_dist = zeros(num_repeats,1);
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
    % Pick alpha indices at random to corrupt
    rand_idx = randperm(num_far_final);
    rand_idx = rand_idx(1:alpha);
    rand_idx_columns = zeros(k,alpha);
    % Corruption of normal nodes
    % Uniform noise
    beta_normal = beta1_min + (beta1_max - beta1_min) * rand(size(F));
    random_sign = sign(rand(size(F)) - 0.5);
    if central_corrupt==1
        F_corrupted = (1 + beta_normal.*random_sign).*F;
    else
        F_corrupted = (1 + beta_normal.*random_sign).*F;
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
        random_sign = sign(rand(size(F1)) - 0.5);
        F_corrupted(rand_idx_1,rand_idx(j)+num_near_final) = (1+beta_corrup.*random_sign).*F1;
    end
    central = anchors(:,m);
    % Re-define anchors removing central node
    anchors = anchors(:,1:m-1);
    % Create the system Aq=b and Aq = b_corrupted
    A = (anchors-central)';
    anchor_norms = vecnorm(anchors).^2;
    % Regularization parameter lambda
    lamda_array = [1 2 5 10 20];
    S_norms = zeros(length(lamda_array),num_near_final+num_far_final);
    % Use error correction to estimate S
    A_1 = [A ones(m-1,1)];
    C = null(A_1')';
    S1 = zeros(m-1,num_far_final+num_near_final);
    target=[target_near target_far];
    % Solve least squares problems 
    N = num_near_final + num_far_final;
    norm_c2 = norm(central)^2;
    if central_corrupt==1
    % Distance from target to central node
        t_c_dist = F_corrupted(m,1:N);
    else
        diffs = target(:,1:N)-central;
        t_c_dist = sum(diffs.^2,1);
    end
    B = 0.5 * ( t_c_dist - F_corrupted(1:m-1,1:N) + anchor_norms' - norm_c2);
    S1 = C\(C*B);
    % Solve the optimization program
    for p = 1:length(lamda_array)
        cvx_begin quiet
        cvx_solver mosek
        cvx_precision high
        variable S(m,num_near_final+num_far_final)
        variable L(m,num_near_final+num_far_final)
        minimize norm_nuc(L)+lamda_array(p)*sum(norms(S)) 
        L+S==F_corrupted;
        cvx_end
        S_norms(p,:) = vecnorm(S,2);
    end
    % Find intersection of indices corresponding to the top k_2
    % values in each row of s_norms
    max_indices = zeros(p,alpha);
    for jj = 1:p
        [~,idx] = maxk(S_norms(jj,:),alpha);
        max_indices(jj,:) = sort(idx);
    end
    S_norms = vecnorm(S1,2);
    [~,idx_s_norms] = maxk(S_norms,alpha);
    identification_error1(i) = 1-(length(intersect(idx_s_norms,rand_idx+num_near_final))/alpha);
    identification_error(i) = 1-(length(intersect(rand_idx+num_near_final,top_repeart(max_indices))')/alpha);
    [~,idx_highest_col_norm] = maxk(vecnorm(F_corrupted),alpha);
    naive_error(i) = 1-(length(intersect(idx_highest_col_norm,rand_idx+num_near_final))/alpha);
    % Apply sparse+dense algorithm to identified nodes to remove  outliers
    rel_err = zeros(alpha,1);
    rel_err2 = zeros(alpha,1);
    anchors_dist = zeros(alpha,1);
    anchors_dist2 = zeros(alpha,1);
    pos_err = zeros(alpha,1);
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
        q_est= robust_alg(A,b) ;
        % Performance metrics
        target_pos= target_far(:,rand_idx(tt));
        pos_err(tt) = norm(target_pos-q_est)^2;
        rel_err(tt) = norm(target_pos-q_est)/norm(target_pos);
        t_a_vec = q_est-[anchors central];
        tar_idx = rand_idx(tt)+num_near_final;
        dist_err(tt) = sum((sqrt(F(1:m,tar_idx))-sqrt(F_corrupted(1:m,tar_idx))).^2);
        anchors_dist(tt) = mean(vecnorm(t_a_vec))/mean(sqrt(F(1:m,tar_idx)));
        % Save for visualization
        if opts.visualize==1
            far_nodes(:,tt) = q_est;
        end
    end
    mean_rel_err(i) = mean(rel_err); 
    mean_pos_err(i) = mean(pos_err); 
    mean_dist_err(i) = sum(dist_err)/(alpha*m);
    mean_anchors_dist(i) = mean(anchors_dist);
    A_unit = A./vecnorm(A);
    % Visualize
    if opts.visualize==1
        anchors_central = [anchors central];
        figure
        scatter(anchors_central(1,:),anchors_central(2,:),'b')
        hold on
        scatter(far_nodes(1,:),far_nodes(2,:),'g','filled')
        hold on
        scatter(target_far(1,rand_idx),target_far(2,rand_idx),'r*')
        labels = {'1','2','3','4'};
        for qq = 1:alpha
            text(far_nodes(1,qq), far_nodes(2,qq), labels{qq}, 'FontSize', 8, 'Color', 'k', 'VerticalAlignment', 'bottom');
        end
        for qq = 1:alpha
            text(target_far(1,rand_idx(qq)), target_far(2,rand_idx(qq)), labels{qq}, 'FontSize', 8, 'Color', 'k', 'VerticalAlignment', 'bottom');
        end
        legend('anchors','estimated positions','true positions')
    end

end
