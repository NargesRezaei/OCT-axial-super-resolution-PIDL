function r = refl (Ni)
% This is a function to compute reflection coeficient of sample layers fot grand trouth plot
    r =  abs(Ni(2:end)-Ni(1:end-1))./(Ni(2:end)+Ni(1:end-1));
end
