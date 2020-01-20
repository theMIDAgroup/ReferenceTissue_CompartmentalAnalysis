% Arows = getrows(A,rows) get the rows of matrix A whose indices are given
% in the vector rows.

function Arows = getrows(A,rows)

Arows=A(rows,:);