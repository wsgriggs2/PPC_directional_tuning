function DRmatC = dataproc_func_cpca(TrainData,TrainLabels,m, ...
    Prior,EvalKeep,DFEmethod)

%DRmatC = dataproc_func_cpca(TrainData,TrainLabels,m, ...
%    Prior,EvalKeep,DFEmethod)
%
%%Input arguments:
%1) TrainData (n x N) matrix where n is the number of samples (instances) and
%N is the dimension of the data. For small sample size problems N >> n.
%2) TrainLabels (n x 1) vector of class labels corresponding to TrainData
%3) m (1 x 1) is the final dimension of the extracted subspaces
%4) Prior (string) specifies the way prior (class) probabilities are
%calculated in the method. Choices: 'empirical' (priors match the,relative
%frequencies of classes), or 'uniform' (priors are uniform-default)
%5) EvalKeep (cell) specifies which eigenvalues/eigenvectors to keep.
%Choices: {'spectrum'} (keep those eigenvalues before the sharpes dip in the
%spectrum), {'median'} (keep eigenvalues higher than the median of nonzero 
%eigenvalues), {'energy',XX} (where 0<XX<=1, keep eigenvalues that capture 
%XX*100% of the cummulative energy--XX to be chosen by the user), 'mean' (keep 
%eigenvalues higher than the median of nonzero eigenvalues-default)
%6) DFEmethod (string) specify which discriminant feature extraction method to
%use. Choices: 'aida' (approximate information discriminant analysis,
%'lda' (linear discriminant analysis), 'ida' (information discriminant analysis) 
%or 'identity' (no further feature extraction, i.e. use identity matrix-default)
%
%Output arguments:
%1) DRmatC (1 x C) cell, whrere C is the number of classes.

%Author: Po T. Wang
%Modified by: Zoran Nenadic DSc, June 18, 2010

Nobs = length(TrainLabels);
if Nobs ~= size(TrainData,1)
    error('Number of observations from TrainData and TrainLabels disagree.');
end

classes = unique(TrainLabels);

Nclass = length(classes);

for c = 1:Nclass
    NtrialA(c) = length(find(TrainLabels == classes(c)));
end

for c = 1:Nclass
    
    idc = find(TrainLabels==classes(c));
    sampmu{c} = mean(TrainData(idc,:),1);
    [coeff, latent] = dataproc_func_princomp(TrainData(idc,:));
    ind = find(latent>0);
    latent = latent(ind);
   
    switch EvalKeep{1}
        case 'median'
            coeffrC{c} = coeff(:,find(latent > median(latent)));
        case 'spectrum'
            [mtemp, indtemp] = max(abs(diff(latent)));
            coeffrC{c} = coeff(:,1:indtemp);
        case 'energy'
            indtemp = find(cumsum(latent) > EvalKeep{2} * sum(latent));
            coeffrC{c} = coeff(:,1:indtemp);
        otherwise
            coeffrC{c} = coeff(:,find(latent > mean(latent)));
    end
    
%     disp([num2str(size(coeffrC{c},2)) ' components kept'])
    
end

sampmuall = mean(TrainData,1);

% Calculate between-class covariance
for c = 1:Nclass
    Data_b(c,:) = sqrt(NtrialA(c) / Nobs) * (sampmu{c} - sampmuall);
end

W_b = dataproc_func_princomp(Data_b);

% Calculate principal subspace basis
switch Prior
    case {'empirical'}
        P_i = NtrialA./Nobs;    %prior class probabilities
    otherwise   %uniform
        P_i = ones(1,Nclass);
        P_i = P_i/sum(P_i);
end

for c = 1:Nclass
    TempBasis = orth([coeffrC{c} W_b]);
    TrainDataProj = TrainData * TempBasis;
    switch DFEmethod
        case {'aida'}
            %download aida.m; sqrtmat.m; invsqrtmat.m; and logmat.m from
            %http://cbmspc.eng.uci.edu/SOFTWARE/AIDA/aida.html
            TDFE = aida(TrainDataProj,TrainLabels,P_i,m);
        case {'lda'}
            TDFE = linear_disc_analysis(TrainDataProj,TrainLabels,m);
        case {'ida'}
            %download ida_feature_extraction_matrix.m; 
            %ida_feature_extraction_matrix.m; orthonormalize.m; and 
            %negative_mu_cg.m from http://cbmspc.eng.uci.edu/SOFTWARE/IDA/ida.html
            %download conjugate_gradient.m; and linesearch.m from:
            %http://www2.imm.dtu.dk/~hbn/Software/conj_grad.m
            %http://www2.imm.dtu.dk/~hbn/Software/linesearch.m
            %Change alpha into Alpha throughout conj_grad.m
            %you can change the parameters of ida directly below
            TDFE = ida_feature_extraction_matrix(m,TrainDataProj, ...
                TrainLabels,'cg',10^(-8)*[1 1],2000,'lda',10);
        otherwise
            TDFE = eye(size(TempBasis,2));
    end
    DRmatC{c} = TempBasis * TDFE';
end
