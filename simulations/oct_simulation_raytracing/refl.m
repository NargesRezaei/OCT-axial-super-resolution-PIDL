function r = refl (Ni)
% This function calculate fresnel reflection coeficient for GT plot
r =  abs(Ni(2:end)-Ni(1:end-1))./(Ni(2:end)+Ni(1:end-1));
end
