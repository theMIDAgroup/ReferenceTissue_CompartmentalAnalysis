function [k1rx,k1x,k2x,k3x,k4x,relerr,nit,Cx,Cxt] = ...
    reconstruction_RTM_realdata(Ct,Ctr_func,Vbt,Vbr,sigma,t,toll,r,ngl,k1r,k1,k2,k3,k4)

% needed variables
alpha = [1,1];
[nodes,weights] = gauss_legendre(ngl); % ngl = 8;
kij = [k1r; k1; k2; k3; k4];

%% initial guess

k = num2cell([k1r,k1,k2,k3,k4]+(rand(1,5)-1/2).*[k1r,k1,k2,k3,k4]);
[k1rx,k1x,k2x,k3x,k4x] = deal(k{:});
% k1rx = rand(1); k1x = rand(1); k2x = rand(1); k3x = rand(1); k4x = rand(1);

x = [k1rx; k1x; k2x; k3x; k4x];

% Solve DIRECT PROBLEM (with initial parameters)
Ax = [-(k2x+k3x) k4x; k3x -k4x];
Cax = concentration_Ca_Ctr(k1rx,Vbr,sigma,Ctr_func,0,0,t,nodes,weights);
Cax = @(tt)(interp1([0 t],[0 Cax],tt,'linear',0));
Cx = concentration_TT(k1x,Ax,Cax,0,[0;0],t,nodes,weights);
Cxt = ( (1 - Vbt)*alpha*Cx + Vbt*Cax(t) ).';

% Relative error
relerr = norm(Cxt-Ct)/norm(Ct);

%% GAUSS-NEWTON Regularized method

% Iteration numbers
nit = 0; nit_max = 30; cont = 0;

while (relerr>toll || abs(k1r-k1rx)/k1r>1.25e-1 || abs(k1-k1x)>2.5e-1) || nit == 0 
    
    % Count iteration's number
    nit = nit+1;
    
    % first entry of D, derivative with respect to k1r
    der_k1r = dF_dk1r(k1rx,Vbr,sigma,Vbt,k1x,Ax,Ctr_func,t,nodes,weights);
    % second entry of D, derivative with respect to k1
    der_k1 = dF_dk1(Vbt,Ax,Cax,t,nodes,weights);
    % derivative with respect to A
    der_A = dF_dA(Vbt,k1x,Ax,Cax,t,nodes,weights);
 
    D = [der_k1r, der_k1, -der_A(:,1), der_A(:,2)-der_A(:,1), der_A(:,3)-der_A(:,4)];    
    g = Ct-Cxt;
        
    % solve
    h = (r*diag([1,1,1,1,1])+D.'*D)\(D.'* g);
    
    % Check for positive coefficients
    xph = x+h;
    if any(xph<=0)
        h(xph<=0)=0;
    end
    x = x+h;
    
    % Refresh the parameters
    k1rx = x(1); 
    k1x = x(2); k2x = x(3); k3x = x(4); k4x = x(5);
    
    % Solve DIRECT PROBLEM --> refresh data
    Ax = [-(k2x+k3x) k4x; k3x -k4x];
    Cax = concentration_Ca_Ctr(k1rx,Vbr,sigma,Ctr_func,0,0,t,nodes,weights);
    Cax = @(tt)(interp1([0 t],[0 Cax],tt,'linear',0));
    Cx = concentration_TT(k1x,Ax,Cax,0,[0;0],t,nodes,weights);
    Cxt = ( (1 - Vbt)*alpha*Cx + Vbt*Cax(t) ).';

    % Relative error
    [relerrprec,relerr] = deal(relerr,norm(Cxt-Ct)/norm(Ct));
    
    %.................................................................%
    if  ( nit>=15 && relerr>=0.3 ) || ( relerr>=0.3 && abs(relerr-relerrprec)<1e-3 ) || ( nit>=25 && abs(relerr-relerrprec)<1e-4 ) || (nit>=nit_max)

        cont = cont+1; nit = 0;
               
        % initial guess
        k = num2cell([k1r,k1,k2,k3,k4]+(rand(1,5)-1/2).*[k1r,k1,k2,k3,k4]);
        [k1rx,k1x,k2x,k3x,k4x] = deal(k{:});
        % k1rx = rand(1); k1x = rand(1); k2x = rand(1); k3x = rand(1); k4x = rand(1);

        x = [k1rx; k1x; k2x; k3x; k4x];
        
        % direct problem
        Ax = [-(k2x+k3x) k4x; k3x -k4x];
        Cax = concentration_Ca_Ctr(k1rx,Vbr,sigma,Ctr_func,0,0,t,nodes,weights);
        Cax = @(tt)(interp1([0 t],[0 Cax],tt,'linear',0));
        Cx = concentration_TT(k1x,Ax,Cax,0,[0;0],t,nodes,weights);
        Cxt = ( (1 - Vbt)*alpha*Cx + Vbt*Cax(t) ).';
        
        % relative error
        relerr = norm(Cxt-Ct)/norm(Ct);
        
        continue
        
    end
    
    %.................................................................%
    if cont>=15
        
        toll = toll+toll*5e-2;
        
        cont = 0; nit = 0;
        
        % Initial guess
        k = num2cell([k1r,k1,k2,k3,k4]+(rand(1,5)-1/2).*[k1r,k1,k2,k3,k4]);
        [k1rx,k1x,k2x,k3x,k4x] = deal(k{:});
        % k1rx = rand(1); k1x = rand(1); k2x = rand(1); k3x = rand(1); k4x = rand(1);

        x = [k1rx; k1x; k2x; k3x; k4x];
        
        % direct problem
        Ax = [-(k2x+k3x) k4x; k3x -k4x];
        Cax = concentration_Ca_Ctr(k1rx,Vbr,sigma,Ctr_func,0,0,t,nodes,weights);
        Cax = @(tt)(interp1([0 t],[0 Cax],tt,'linear',0));
        Cx = concentration_TT(k1x,Ax,Cax,0,[0;0],t,nodes,weights);
        Cxt = ( (1 - Vbt)*alpha*Cx + Vbt*Cax(t) ).';

        % relative error
        relerr = norm(Cxt-Ct)/norm(Ct);
        
        continue
        
    end
    
end