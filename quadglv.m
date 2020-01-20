% intf = quadglv(f,a,b,x,w) approximates the integral of the complex 
% array-valued function f from a to b using the Gauss-Legendre quadrature 
% rule of order length(x)=length(w).
%
% f is a function handle. The function A=fun(t) should accept a scalar
% argument t and return an array result A, whose components are the
% integrands evaluated at t.
%
% x and w are respectively the nodes and weights of the Gauss-Legendre
% quadrature on [-1,1] given by the function gauss_legendre.

function intf = quadglv(f,a,b,x,w)

% Nodes and weigths of the Gauss-Legendre quadrature rule on [a,b]
mid=(b+a)/2;
hl=(b-a)/2;
% Nodes
t=hl*x+mid;
% Weights
wt=hl*w;

% f at the nodes t, result arrays in cells
intf=cellfun(f,num2cell(t),'UniformOutput',false);

% Multiply each array by the corresponding weight
intf=cellfun(@times,num2cell(wt),intf,'UniformOutput',false);

% Concatenate arrays along a further dimension
nf=ndims(intf{1})+1;
intf=cat(nf,intf{:});

% Sum along the added dimension to have the approximation of the integral
intf=sum(intf,nf);

end