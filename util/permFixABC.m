function [M_perm, N_perm, P_perm] = permFixABC(M, N, P, r)
% Generate partially permuted data triples [A1_perm, B1_perm, C1_perm]
% Input:
%       M in 4 x 4 x 1
%       N in 4 x 4 x n
%       P in 4 x 4 x n
%       r : scrambling rate
% Output:
%       M_perm in 4 x 4 x n (n copies of M)
%       N_perm in 4 x 4 x n (permuated N)
%       P_perm in 4 x 4 x n (samme as P)


%%
n = size(N, 3);

M_perm = zeros(4, 4, n);
for i = 1:n
    M_perm(:,:,i) = M;
end

N_perm = scrambleData(N, r);

P_perm = P;

end
