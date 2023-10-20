function w = Lambert_W(branch, x)
% Lambert_W Lambert's W function.
%    W = lambertw(Z) solves W*exp(W) = Z.
%    W = lambertw(K,Z) is the K-th branch of this multi-valued function.
 
%    References:
%    [1] Robert M. Corless, G. H. Gonnet, D. E. G. Hare,
%    D. J. Jeffrey, and D. E. Knuth, "On the Lambert W Function",
%    Advances in Computational Mathematics, volume 5, 1996, pp. 329-359.
 
%    [2] Corless, Jeffrey, Knuth, "A Sequence of Series
%    for The Lambert W Function", ISSAC '97

%    Copyright Lateef Adewale Kareem 2022.


if nargin < 2
    x = branch;  branch = 0;
end
if numel(x) > 1 && numel( branch) == 1
    branch = repmat(branch, size(x));
end
if numel(x) == 1 && numel( branch) > 1
    x = repmat(x, size(branch));
end

if numel(x) > 1 && numel( branch) > 1
    sx = size(x); sb = size(branch);
    if(~isequal(sx, sb))
        if(sx(1) == 1 && sb(2) == 1)
            x = repmat(x, sb(1), 1);
            branch = repmat(branch, 1, sx(2));
        elseif (sx(2) == 1 && sb(1) == 1)
            branch = repmat(branch, sx(1), 1);
            x = repmat(x, 1, sb(2));
        end
    end
end




% Effective starting guess
w = ones(size(x));
w(branch ~= 0) = -2;

% lambert w for branches 0 and -1
indx = branch == 0 | branch == -1;
w(indx) = Haley(w(indx), x(indx), 0);

% lambert w for other branches
indx = ~(branch == 0 & branch ==-1);
w(indx) = Haley(w(indx), x(indx), branch(indx));

% correct for x = 0 and -1/e
w(x==0) = 0;
w(x==-1/exp(1)) = -1;
end

function w = Haley(w, x, branch)
    % Haley's method
    IsConverged = 0;
    while (~IsConverged)
       v = w;
       f = w + log(w) - log(x) - 2*pi*1i*branch;
       fp = 1 + 1./w;
       fpp = -1./(w.*w);
       % Iterate to make this quantity zero
       w = w - f ./ (fp - f .* fpp ./ (2 * fp));
       % check convergence
       IsConverged = all(abs(w - v)./abs(w) < 1.e-15);
    end
end