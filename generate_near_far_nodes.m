function [anchors, target_near, target_far] = generate_near_far_nodes(r,m,opts)
    % 
    num_near = opts.num_near;
    num_far = opts.num_far;    
    % Generate random points in the interval [-400,400]. Apply K-means
    % with m clusters, and designate the centers as anchors
    P = -400+800.*rand(r,1000);
    % K-means to initialize anchors
    [~,Centers]=kmeans(P',m);
    anchors = Centers';
    % Set up near target nodes
    target_near = -400+800*rand(2,num_near);
    % Reject target nodes that are too close to each other
    % Reject target nodes that are too close to the anchors
    nodes_near = [target_near anchors];
    near_target_anchors_distances = squareform(pdist(nodes_near'));
    near_target_anchors_distances = near_target_anchors_distances(1:num_near,1:end);
    flag_close = near_target_anchors_distances<40;
    idx_far = sum(flag_close,2)==1;
    target_near = target_near(:,idx_far);
    % construct far targets
    target_far1x = -1200+200*rand(1,num_far);
    target_far1y = -600+1000*rand(1,num_far);
    target_far1 =[target_far1x;target_far1y];
    target_far2x = 1200+200*rand(1,num_far);
    target_far2y = -600+1000*rand(1,num_far);
    target_far2 =[target_far2x;target_far2y];
    target_far = [target_far1 target_far2];
    % Reject target nodes that are too close to each other
    % Reject target nodes that are too close to the anchors
    nodes_far = [target_far anchors];
    far_target_anchors_distances = squareform(pdist(nodes_far'));
    far_target_anchors_distances = far_target_anchors_distances(1:2*num_far,1:end);
    flag_close = far_target_anchors_distances<60;
    idx_far = sum(flag_close,2)==1;
    target_far = target_far(:,idx_far);
end






