% LOGAN PLOT on Ca and Ctr from t0!=0: 
% int_t0^t Ctr / Ctr = alpha * int_t0^t Ca / Ctr + beta
% alpha = (1-Vbr) * k1r/k2r + Vbr
% where - x_logan = int_t0^t Ca / Ctr
%       - y_logan = int_t0^t Ctr / Ctr
%       - t0 is the starting point for the linear regression
%       - t starts from t1
%       - p_logan = [slope,intercept]
    
function [slope,x_logan,y_logan] = logan_plot_Ctr(Ca,Ctr,t0,t,x,w)
    
    nt = length(t);

    Q_Ctr = zeros(1,nt);
    Q_Ca = zeros(1,nt);
    
    f_Ctr = @(u)( Ctr(u) );
    Q_Ctr(1) = quadglv(f_Ctr,t0,t(1),x,w); % if t(1) = t0, Q = 0 (int_t0^t0 = 0 OK!)

    f_Ca = @(u)( Ca(u) );
    Q_Ca(1) = quadglv(f_Ca,t0,t(1),x,w); % if t(1) = t0, Q = 0 (int_t0^t0 = 0 OK!)

    for n=2:nt
        Q_Ctr(n) =  Q_Ctr(n-1) + quadglv(f_Ctr,t(n-1),t(n),x,w);
        Q_Ca(n) =  Q_Ca(n-1) + quadglv(f_Ca,t(n-1),t(n),x,w);
    end
      
    x_logan = Q_Ca./Ctr(t);
    y_logan = Q_Ctr./Ctr(t);
        
    p_logan = polyfit(x_logan,y_logan,1);
    slope = p_logan(1);

end