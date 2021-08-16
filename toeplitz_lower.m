
function B = toeplitz_lower(A,nr,nc,m)

%+----------------------------------------------------------------+%
%                                                                  %
% Generate a lower triangular Toeplitz matrix from a block column  %
%          matrix, containing the elements of this Toeplitz matrix %
% Inputs: A - block column matrix, containing the elements of this %
%             Toeplitz matrix, size(A)=m*nr x mc                   %
%         nr,nc - the row and column dimension of an individual    %
%                 block element                                    %
%         m - number of block elements in A and B                  %
% Outputs:B - the lower triangular Toeplitz matrix                 %
%                                                                  % 
% Author: Jianfei Dong (c), Jan 7, 2008                            %
% Delft Center for Systems and Control, Delft Univ of Technology   %
%                                                                  %
%+----------------------------------------------------------------+%


nin = nargin;
if nin < 4
    error('Wrong number of input arguments');
end

[r,c] = size(A);
if c ~= nc || r < m*nr
    error('Wrong numbers of rows or columns. A must be a block column vector.');
end

B = []; %zeros(m*nr,m*nc);
for k = 1:m
    B = [ B, [ zeros((k-1)*nr,nc); A(1:end-(k-1)*nr,:) ] ];
end
