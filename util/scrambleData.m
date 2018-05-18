function M_perm = scrambleData(M, s_rate)
% Permutation function
% Input:
%       M: set of SE3 matrices
%       s_rate: scrambling rate

M_index = (1:size(M, 3));

for i = 1:length(M_index)
    
    if rand <= 0.01*s_rate
        index = randi(size(M, 3), 1);
        M_index([i index]) = M_index([index i]);
    end
    
end

M_perm = M(:, :, M_index);

end