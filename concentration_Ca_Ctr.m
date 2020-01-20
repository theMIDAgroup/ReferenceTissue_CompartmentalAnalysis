% Ca = concentration_Ca_Ctr(k1r,Vbr,sigma,Ctr,t0,Ca0,t,x,w) computes the 
% Ca concentration starting from 
% Ctr = (1-Vbr)*Cr+Vbr*Ca with the kinetic parameter k1r and 
% the adimensional parameter sigma(=(1-Vbr)/Vbr+k2r/k1r)

% concentration_Ca_Ctr allow directly the computation of the concentration 
% Cr starting from t0 != 0

% Ca = 1/Vbr * [ Ctr - (1-Vbr)/Vbr * k1r * int_t0^t ( exp(-k1r*sigma*(t-u)) * Ctr(u) ) du ] + exp(-k1r*sigma*(t-t0)) * Ca0

function Ca = concentration_Ca_Ctr(k1r,Vbr,sigma,Ctr,t0,Ca0,t,x,w)

    nt = length(t);
    vr = (1-Vbr)/Vbr;    
    Ca_int = zeros(1,nt);
    
    if isa(Ctr,'function_handle') == 0
        Ctr = @(tt)(interp1([0 t],[0 Ctr],tt,'linear',0));
    end
    
    f = @(u)( 1/Vbr * vr*k1r * Ctr(u) * exp(-k1r*sigma*(t(1)-u)) );
    Ca_int(1) = exp(-k1r*sigma*(t(1)-t0)) * Ca0 + quadglv(f,t0,t(1),x,w);

    for n=2:nt 
        
        f = @(u)( 1/Vbr * vr*k1r * Ctr(u) * exp(-k1r*sigma*(t(n)-u)) );
        Ca_int(n) = exp(-k1r*sigma*(t(n)-t(n-1))) * Ca_int(n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
    Ca = 1/Vbr * Ctr(t) - Ca_int;
    
end