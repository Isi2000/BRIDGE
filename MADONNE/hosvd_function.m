function [hatT,U,S,sv,nn]=hosvd_function(Tensor,varepsilon1,nn,n,TimePos)

% PERFORM HOSVD AND RETAINS ALL THE SINGULAR VALUES
% [L,I,J,K]=size(Tensor)
% n=[L,I,J,K];
[TT S U sv n] = hosvd(Tensor, n);

% SET THE TRUNCATION OF THE SINGULAR VALUES USING varepsilon1 (automatic
% truncation)
%nn=[L,0,0,0];

for i=2:length(nn)
    count=0;
    for j=1:size(sv{i},1)
        if sv{i}(j)/sv{i}(1)>=varepsilon1
            count=count+1;
        else
            break
        end
    end
    nn(i)=count;
end

('Initial number of singular values')
n
('Number of singular values retained')
nn

% HOSVD retaining n singular values: reconstruction of the mdoes
[TT S U sv nn] = hosvd(Tensor, nn);

% Construct the reduced matrix containing the temporal modes
for pp=1:length(nn)
    UT{pp}=U{pp}';
end

for kk=1:nn(TimePos)
    hatT(kk,:)=sv{TimePos}(kk)*UT{TimePos}(kk,:);
end