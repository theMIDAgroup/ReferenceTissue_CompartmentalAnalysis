% Acols = getcols(A,cols) get the columns of matrix A whose indices are given
% in the vector cols.

function Acols = getcols(A,cols)

Acols=A(:,cols);
