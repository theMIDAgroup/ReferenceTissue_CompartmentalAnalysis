function ind_t0 = compute_logan_t0(t,Ca,Ctr)

    aux = Ca(t)./Ctr(t);
    aux_diff = abs(diff(aux));

    ind_t0 = find(abs(aux_diff)<3*10^(-2),1);

    if isempty(ind_t0)    
        ind_t0 = length(t)-floor(length(t)/6);
        % error('WARNING: lambda not found - hand checking')    
    end

end