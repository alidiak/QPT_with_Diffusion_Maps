function [centroids, idx] = k_means_func(X, K, max_iters) %, vargin)
% kmeans runs the K-Means algorithm on data matrix X, where each row of X
% is a single sample
% [centroids, idx] = kmeans(X,K,max_iters)  will automatically initiallize
% clusters based on the domain, range of the data X or will initialize them
% by randomly chosing one of the data points to start at.The kmeans
% algorithm will then run for n=max_iters and with number of clusters equal
% to K. The algorithm will return data space position of clusters in
% centroids as well as the indices assigned to the data based on their
% nearest clusters.
%

[m n] = size(X); % data m is sample size, n is data dimension

%% Initialize the centroids

% Randomly reorder the indices of examples
randidx = randperm(size(X, 1));

% Take the first K examples as centroids
centroids = X(randidx(1:K), :);

%% Run the k-means algorithm
%  previous_centroids = centroids;
idx = zeros(m, 1);

% Run K-Means
for i=1:max_iters
    
    % Output progress
    fprintf('K-Means iteration %d/%d...\n', i, max_iters);
    
    % For each sample in X, assign it to the closest centroid
    distance = zeros(K,1);
    for i=1:m
    
        for k = 1:K
            % calculating euclidian distance squared between x(i) and u_k
            % for each sample i, calculate the distance to each centroid k
            % then pick the index idx to be the centroid it is closest to
            distance(k) = sum((X(i,:)-centroids(k,:)).^2);   
        end
        % find the index of the centroid that is least distant and assign x(i)
        % to it
        [M,idx(i)] = min(distance);
    end
    
    % Optionally, plot progress here
%     if plot_progress
%         plotProgresskMeans(X, centroids, previous_centroids, idx, K, i);
%         previous_centroids = centroids;
%         fprintf('Press enter to continue.\n');
%         pause;
%     end
    
    % Given the memberships, compute new centroids
    for k= 1:K
        % idx == k returns a binary array whereever the index of the closest 
        % centroid matches k. Inputing this into X(,:) then returns the
        % coordinates of those sites closest to centroid k.
        points_nearest_k = X(idx==k,:);
   
        % compute the new centroid position by taking the mean of the
        % samples assigned to it
        centroids(k,:)= mean(points_nearest_k);    
    end
end

% Hold off if we are plotting progress
% if plot_progress
%     hold off;
% end

% Plot final results
if n==2
    figure
    palette = hsv(K + 1);   colors = palette(idx, :); % assign K colors
    scatter(X(:,1), X(:,2), 15, colors); hold on
    plot(centroids(:,1), centroids(:,2), 'x', ... % plot final centroids
    'MarkerEdgeColor','k', ...
    'MarkerSize', 10, 'LineWidth', 3);
end

end


%% Testing stuff

% calculating the average within cluster standard deviation
% avg_cluster_dev=0; avg_intercluster_dist=0;
% for k=1:K
%     % X(idx==k,:) returns the samples belonging to cluster k, then, the 
%     % sqrt of the average of the distance squared is the standard deviation.
%     cluster_dev = sqrt(mean((X(idx==k,:)-centroids(k,:)).^2,1));
%     avg_cluster_dev=avg_cluster_dev+cluster_dev;
%     
%     % Calculating the average intercluster distance
%     for k2=1:K
%         D=sqrt(sum((centroids(k,:)-centroids(k2,:)).^2));
%         k
%         k2
%         D
%         avg_intercluster_dist=avg_intercluster_dist+D;
%     end
% end
% avg_cluster_dev=avg_cluster_dev/K;
% % possibly 0.5 times this as distances are counted twice
% avg_intercluster_dist=(1/(K*(K-1)))*avg_intercluster_dist; 
