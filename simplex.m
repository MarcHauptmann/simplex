function x = simplex(c, A, b)
%SIMPLEX Simplex method
%   x = SIMPLEX(c, A, b)
%
%   solves the problem
%   $$ \min_x c^Tx $$
%   s.t. Ax = b
%        x >= 0
%

    [m, n] = size(A);
    BIndices = 1:m;
    
    while 1
        NIndices = setdiff(1:n, BIndices);

        B = A(:,BIndices);
        N = A(:,NIndices);

        x = zeros(n, 1);
        s = zeros(n,1);

        x(BIndices) = B\b;

        lambda = B'\c(BIndices);
        s(NIndices) = c(NIndices) - N'*lambda;

        if all(s>=0)
            fprintf('Lösung gefunden');
            return;
        end

        [sq, q] = min(s);

        d = B\A(:,q);

        [xqp, p] = min(x(BIndices)./d);

        x(BIndices) = x(BIndices) - d*xqp;
        x(q) = xqp;
        
        BIndices = [setdiff(BIndices, BIndices(p)), q];
    end
end
