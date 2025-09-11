function [TT S U sv n] = hosvd(T, eps)
%HOSVD High Order SVD of a multidimensional array
%	[TT S U sv tol] = HOSVD(T)
%	[TT S U sv tol] = HOSVD(T, n)
%
%
%	T    - multidimensional array
%	n    - SV 1..n(i) are taken. If not
%
%	S    - decomposed core tensor so that T==tprod(S, U)
%	U    - matrices for each dimension
%	sv   - n-mode singular values
%	tol  - largest dropped singular value (hosvd truncates small sv)
%
%	eg. [TT, S, U, sv, tol] = hosvd(ones(3,4,5), [2 2 3])
%
%	See also TPROD, SVDTRUNC

M = size(T);
P = length(M);

U = cell(1,P);
UT = cell(1,P);
sv = cell(1,P);
producto=numel(T);

for i = 1:P
    n(i)=min(M(i),producto/M(i));
    A = ndim_unfold(T, i);

    % SVD based reduction of the current dimension (i)
    [Ui svi] = svdtrunc(A, n(i));
    
    if n(i)<2
        U{i}=[Ui(:,1),zeros(size(Ui(:,1)))];
    else
    count=0;
    for j=1:size(svi,1)
        if svi(j)/svi(1)>=eps
            count=count+1;
        else
            break
        end
    end
    n(i)=count;          
    U{i} = Ui(:,1:n(i));
    end
    %U{i} = Ui;
    UT{i} = U{i}';
    sv{i} = svi(1:n(i));
    %tol(i) = toli;
end

%   vega=size(UT{2})
%   vega=size(UT{3})
S = tprod(T, UT);
TT=tprod(S, U);

