% [lambda_GCV,V_lambda] = GCV(A,g,vect_lambda)
% The Generalized Cross Validation (GCV) computes the function V_lambda in
% correspondence of the values vect_lambda. V_lambda depends only on the
% SVD of the matrix of the problem A (zero computational cost).
%
% A is the matrix of the problem  | ==> Af=g
% g is the (noisy) data           |
% [HP: A is overdetermined, rows>columns]
% vect_lambda is the vector of the possible values for lambda (--> the regularization parameter for Tikhonov)
% ==>  Tychonov regularization: (lambda*I+A'*A)h=A'*g
%
% lambda_GCV is the value of the regularization parameter at
% which V_lambda assumes its minimum.

function [lambda_GCV,V_lambda] = GCV(A,g,vect_lambda)

    V_lambda = zeros(size(vect_lambda));

    m = size(A,1); n = size(A,2);
    [U,S] = svd(A);
    vect_sigma = diag(S);
    hh = U'*(g);

    for l=1:length(vect_lambda)

        V_num=0;
        V_den=0;

        for i=1:n

            V_num = V_num + ( vect_lambda(l)*hh(i) / (vect_lambda(l)+vect_sigma(i)^2) )^2;
            V_den = V_den + ( vect_sigma(i)^2 / (vect_sigma(i)^2 + vect_lambda(l)) );

        end

        V_num = V_num + sum(hh(n+1:m).^2);
        V_lambda(l) = 1/m * ( V_num / (1 - (1/m)*V_den)^2 );

    end

    lambda_GCV = vect_lambda(V_lambda==min(V_lambda)); % lambda value corresponding to the minimum of the fuction V_lambda

    % Control on the uniqness of the minimum: in general, the minimum of
    % the function V_lambda may not be isolated, then lambda_GCV may be a
    % vector of points (with consecutive values).
    % We just pick an entry of this vector, for example the last one 
    % (which is also the bigger one).
    
    if length(lambda_GCV) > 1

        lambda_GCV = max(lambda_GCV);

    end

end