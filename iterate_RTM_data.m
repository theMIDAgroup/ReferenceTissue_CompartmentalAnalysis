function [K1rx,K1x,K2x,K3x,K4x,Cax,Cxt,relerr_TT,relerr_IF,nit,crit,relerrp] = ...
    iterate_RTM_data(glnodes,glweights,t,alpha,Vbt,Vbr,sigma_logan,Ct,Ctr,Ca,...
            K1rx,K1x,K2x,K3x,K4x,Cax,Cxt,relerr_TT,nit,crit)
        
% crit is a 5x5 matrix, each row for the relative errors between successive
% iterations, ie for example, if nit=8:
% first row: relative errors between step 3 and step 4
% second row: relative errors between step 4 and step 5
% third row: relative errors between step 5 and step 6
% fourth row: relative errors between step 6 and step 7
% fifth row: relative errors between step 7 and step 8
% then, after the function, nit=9 and crit becomes
% first row: relative errors between step 4 and step 5
% second row: relative errors between step 5 and step 6
% third row: relative errors between step 6 and step 7
% fourth row: relative errors between step 7 and step 8
% fifth row: relative errors between step 8 and step 9

nit = nit+1;
Ax = [-(K2x+K3x) K4x; K3x -K4x];
x = [K1rx;K1x;K2x;K3x;K4x];

% first entry of D, derivative with respect to k1r
der_k1r = dF_dk1r(K1rx,Vbr,sigma_logan,Vbt,K1x,Ax,Ctr,t,glnodes,glweights);
% second entry of D, derivative with respect to k1
der_k1 = dF_dk1(Vbt,Ax,Cax,t,glnodes,glweights);
% derivative with respect to A
der_A = dF_dA(Vbt,K1x,Ax,Cax,t,glnodes,glweights);

D = [der_k1r, der_k1, -der_A(:,1), der_A(:,2)-der_A(:,1), der_A(:,3)-der_A(:,4)];    
g = Ct-Cxt;

% regularization parameter
% GCV
% vect_r = 1e5:1e5:1e7; vect_r = vect_r';
% [r,~] = GCV(D,g,vect_r);
r=1e7;

% solve
h = (r*diag([1,1,1,1,1])+D.'*D)\(D.'* g);

% Check for positive coefficients
xph = x+h;
if any(xph<=0)
    h(xph<=0)=0;
end
x = x+h;

% keep trace of previous values
[K1rxp,K1xp,K2xp,K3xp,K4xp,relerrp]=deal(K1rx,K1x,K2x,K3x,K4x,relerr_TT);

K1rx=x(1);K1x=x(2);K2x=x(3);K3x=x(4);K4x=x(5);
Ax = [-(K2x+K3x) K4x; K3x -K4x];

Cax = concentration_Ca_Ctr(K1rx,Vbr,sigma_logan,Ctr,0,0,t,glnodes,glweights);
Cax = @(tt)(interp1([0 t],[0 Cax],tt,'linear',0));
Cx = concentration_TT(K1x,Ax,Cax,0,[0;0],t,glnodes,glweights);
Cxt = ( (1 - Vbt)*alpha*Cx + Vbt*Cax(t) ).';

relerr_TT = norm(Cxt-Ct)/norm(Ct);
relerr_IF = norm(Cax(t)-Ca(t))/norm(Ca(t));

crit(1,:)=[];
crit(end+1,:)=[abs(K1rx-K1rxp)/abs(K1rx),...
    abs(K1x-K1xp)/abs(K1x),...
    abs(K2x-K2xp)/abs(K2x),...
    abs(K3x-K3xp)/abs(K3x),...
    abs(K4x-K4xp)/abs(K4x),...
    abs(relerr_TT-relerrp)/abs(relerr_TT)];

end