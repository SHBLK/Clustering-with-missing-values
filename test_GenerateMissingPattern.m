function Miss_matrix = test_GenerateMissingPattern(n,p,missing_prob)

%This function generates the missing value patterns
Miss_matrix=binornd(ones(n,p),(1-missing_prob)*ones(n,p));
while any(sum(Miss_matrix,1)==0) || any(sum(Miss_matrix,2)==0)
    Miss_matrix=binornd(ones(n,p),(1-missing_prob)*ones(n,p));
end

end

