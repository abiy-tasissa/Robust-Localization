function q_estimated_alg= robust_alg(A,b_corrupted)  
    m_2 = length(b_corrupted);
    sz_A = size(A);
    A_tilde= [A ones(sz_A(1),1)];
    RR = null(A_tilde')';
    cvx_begin quiet
    cvx_solver mosek
    cvx_precision high
    variable s(m_2,1)
    minimize sum(abs(s))
    RR*s==RR*b_corrupted;
    cvx_end
    qc = A_tilde\(b_corrupted-s);
    q_estimated_alg = qc(1:end-1);
end