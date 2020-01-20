% C = concentration_TT(k1,A,Ca,t0,C0,t,x,w) computes the concentration for 
% the TT tissue compartments C=[Cf;Cp] with the kinetic parameters 
% k1,k2,k3,k4 solution of
% C'= A*C + k1*Ca*e1 with e1=[1;0] and initial condition C(t0)=C0.

% concentration_TT allow directly the computation of the concentration C
% starting from t0 != 0

% C = k1 * int_t0^t ( exp(A*(t-u)) * Ca(u) * e1 ) du + exp(A*(t-t0)) * C0

function C = concentration_TT(k1,A,Ca,t0,C0,t,x,w)

    nt = length(t);
    C = zeros(2,nt);
    
    if isa(Ca,'function_handle') == 0
        Ca = @(tt)(interp1([0 t],[0 Ca],tt,'linear',0));
    end
        
    f = @(u)( k1 * Ca(u) * getcols(expm((t(1)-u)*A),1) );
    C(:,1) = expm((t(1)-t0)*A) * C0(:) + quadglv(f,t0,t(1),x,w);

    for n=2:nt
               
        f = @(u)( k1 * Ca(u) * getcols(expm((t(n)-u)*A),1) );           
        C(:,n) = expm((t(n)-t(n-1))*A) * C(:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
end