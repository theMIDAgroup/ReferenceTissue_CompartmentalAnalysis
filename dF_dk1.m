function der_k1 = dF_dk1(Vbt,A,Ca,t,x,w)

    alpha = [1-Vbt,1-Vbt];
    
    nt = length(t);
    C = zeros(2,nt);
        
    if isa(Ca,'function_handle') == 0
        Ca = @(tt)(interp1([0 t],[0 Ca],tt,'linear',0));
    end
    
    f = @(u)( Ca(u) * getcols(expm((t(1)-u)*A),1) );    
    C(:,1) = quadglv(f,0,t(1),x,w);

    for n=2:nt
        
        f = @(u)( Ca(u) * getcols(expm((t(n)-u)*A),1) );        
        C(:,n) = expm((t(n)-t(n-1))*A) * C(:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
    
    der_k1 = (alpha * C).';
    
end

