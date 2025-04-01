function Anew = aida(Data,Member,P_i,m)

%Data : [Ntr x n] matrix of data (Ntr - # of trials, n - # of
%attributes)
%Member : [Ntr x 1] matrix of labels
%make sure Data and Member is train data and its respective
%class memberships. 
%p : class probabilities (from training data)
%m : the number of features

no_class = length(P_i); %# of classes
Ntr = size(Data,1);  %total # of trials 
no_feature = size(Data,2);  %total # of feature

M = mean(Data,1);  %total mean
Si=[];PhiB=[];
Sw = 0;
for i = 1:no_class
  ind = find(Member == i-1);
  Mi = mean(Data(ind,:),1); %class-conditional means
  mu_i{i}=Mi';
  PhiB = [PhiB; sqrt(P_i(i))*[Mi - M]]; %note: PhiB [c x n] matrix 
  Si{i}=cov(Data(ind,:),1);
  Sw = Sw + P_i(i) * cov(Data(ind,:),1);
end


    Sb=PhiB'*PhiB;
    St=Sw+Sb;
    Ss = sqrtmat(Sw);
    Sinv = inv(Sw);
    Sis = invsqrtmat(St);
    Sc=0;


   %sphering
for i = 1:no_class

    Zi{i} = Sis * Si{i} * Sis;  %transformed si's

end


%---------------------- NEW Matrix ---------% 

Znew = zeros(no_feature,no_feature);

for i = 1:no_class

  Znew = Znew - P_i(i) * logmat(Zi{i});
 
end


[V, D] = eig(Znew);

[D iD] = sort(diag(D));

V = fliplr(V(:,iD));
V = V(:,1:m);

Anew = V'*Sis;
Anew = Anew/norm(Anew);   %this is the feature extraction matrix




