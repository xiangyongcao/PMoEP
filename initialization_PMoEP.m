function R = initialization_PMoEP(X, k, type)
N = length(X);
switch type
    case 'kmeans'
        [Ind, ~] = kmeans(X,k,'emptyaction','singleton','replicates',1);
        R = full(sparse(1:N,Ind,1,N,k,N));
    case 'random'
        idx = randsample(N,k);%从1到N中随机抽取k个数
        tmp = X(:,idx);
        [~,label] = max(bsxfun(@minus,tmp'*X,dot(tmp,tmp,1)'/2),[],1);%列排序，label是最值得位置
        [u,~,label] = unique(label);
        while k ~= length(u)
            idx = randsample(N,k);
            tmp = X(:,idx);
            [~,label] = max(bsxfun(@minus,tmp'*X,dot(tmp,tmp,1)'/2),[],1);
            [u,~,label] = unique(label);
        end
        R = full(sparse(1:N,label,1,N,k,N));
end