function der_k1r = dF_dk1r(k1r,Vbr,sigma,Vbt,k1,A,Ctr,t,x,w) 
    
    nt = length(t);
    vr = (1-Vbr)/Vbr;
    alpha = [1-Vbt,1-Vbt];
            
    if isa(Ctr,'function_handle') == 0
        Ctr = @(tt)(interp1([0 t],[0 Ctr],tt,'linear',0));
    end
    
    % Cr~ = int_0^t exp(-k1r*sigma*(t-u)) Ctr(u) du
    % Cr* = -(1-Vbr)/Vbr * [Cr~ + k1r * dCr~/dk1r] = 
    %     = -(1-Vbr)/Vbr * int_0^t (1-k1r*sigma*(t-u)) * exp(-k1r*sigma*(t-u)) Ctr(u) du
    Cr_star = zeros(1,nt); 

    f = @(u)( -vr * (1-k1r*sigma*(t(1)-u)) * Ctr(u) * exp(-k1r*sigma*(t(1)-u)) );
    Cr_star(1) = quadglv(f,0,t(1),x,w);

    for n=2:nt         
        f = @(u)( -vr * (1-k1r*sigma*(t(n)-u)) * Ctr(u) * exp(-k1r*sigma*(t(n)-u)) );
        Cr_star(n) = exp(-k1r*sigma*(t(n)-t(n-1))) * Cr_star(n-1) + quadglv(f,t(n-1),t(n),x,w);        
    end
    Cr_star = @(tt)(interp1([0 t],[0 Cr_star],tt,'linear',0));
    
    % dF/dk1r = (1-Vbt)*alpha * k1/Vbr * [int_0^t exp(A*(t-u)) * e1 * Cr*(u) du] +
    %           + Vbt/Vbr * Cr*(t)
    der = zeros(2,nt);
        
    g = @(u)( Cr_star(u) * getcols(expm((t(1)-u)*A),1) );    
    der(:,1) = quadglv(g,0,t(1),x,w);

    for n=2:nt;        
        g = @(u)( Cr_star(u) * getcols(expm((t(n)-u)*A),1) );       
        der(:,n) = expm((t(n)-t(n-1))*A) * der(:,n-1) + quadglv(g,t(n-1),t(n),x,w);        
    end
    
    der_k1r = ( alpha * k1/Vbr * der + Vbt/Vbr * Cr_star(t) ).';
    
end