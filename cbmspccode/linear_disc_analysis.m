function [Tlda] = linear_disc_analysis(TrainData,TrainLabels,m)

[Nt, n] = size(TrainData);
classes = unique(TrainLabels);
Nclass = length(classes);

Mean = zeros(n,1);
Sigma_w = zeros(n,n);
Sigma_b = zeros(n,n);

for i = 1:Nclass
    
    ind = find(TrainLabels == classes(i));
    
    prob_i(i) = length(ind)/Nt;
    Mean_i(:,i) = mean(TrainData(ind,:),1)';
    Sigma_i(:,:,i) = cov(TrainData(ind,:));
    
    %overall mean
    Mean = Mean + prob_i(i) * Mean_i(:,i);
    
    %within class matrix
    Sigma_w = Sigma_w + prob_i(i) * Sigma_i(:,:,i);

    %between class matrix
    Sigma_b = Sigma_b + prob_i(i) * [Mean_i(:,i) * Mean_i(:,i)'];

end

Sigma_b = Sigma_b - Mean * Mean';
Sigma = Sigma_w + Sigma_b;

Rsb = rank(Sigma_b);

mnew = min([Rsb,Nclass-1,m]);

if mnew < m
    disp('warning: your choice of subspace dimension is too high')
    disp('warning: switching to default')
end

[V D] = eig(inv(Sigma) * Sigma_b);
D = diag(D);
[D Ind] = sort(D);
V = V(:,Ind);
V = fliplr(V);
Tlda = V(:,1:mnew)';

