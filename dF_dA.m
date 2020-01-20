function der_A = dF_dA(Vbt,k1,A,Ca,t,x,w)
 
    nt = length(t);
    alpha = [1-Vbt,1-Vbt];
    der = zeros(2,4,nt);
        
    C = concentration_TT(k1,A,Ca,0,[0;0],t,x,w);
    
    Cu = @(u)( concentration_TT(k1,A,Ca,0,[0;0],u,x,w) );    
    f = @(u)( kron((Cu(u)).',expm((t(1)-u)*A)) );   
    der(:,:,1) = quadglv(f,0,t(1),x,w);
    
    for n=2:nt;

        Cu = @(u)( concentration_TT(k1,A,Ca,t(n-1),C(:,n-1),u,x,w) );
        f = @(u)( kron((Cu(u)).',expm((t(n)-u)*A)) );        
        der(:,:,n) = expm((t(n)-t(n-1))*A) * der(:,:,n-1) + quadglv(f,t(n-1),t(n),x,w);
        
    end
   
    der_A = reshape(alpha*reshape(der,2,4*nt),4,nt).';
     
end    
