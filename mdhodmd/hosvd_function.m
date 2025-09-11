function [hatT,U,S,sv,nn,N]=hosvd_function(Tensor,varepsilon1,TimePos)

% PERFORM HOSVD AND RETAINS ALL THE SINGULAR VALUES
n=size(Tensor)
[TT S U sv nn] = hosvd(Tensor,varepsilon1);

% SET THE TRUNCATION OF THE SINGULAR VALUES USING varepsilon1 (automatic
% truncation)

('Initial number of singular values')
n
('Number of singular values retained')
nn

% HOSVD retaining n singular values: reconstruction of the mdoes

% Construct the reduced matrix containing the temporal modes
for pp=1:ndims(Tensor)
    UT{pp}=U{pp}';
end

for kk=1:nn(TimePos)
    hatT(kk,:)=sv{TimePos}(kk)*UT{TimePos}(kk,:);
end

[N,~]=size(hatT);