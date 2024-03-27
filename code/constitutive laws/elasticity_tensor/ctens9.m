%--------------------------------------------------------------------------
% Evaluates the constitutive tensor (in Voigt notation) for material type 8.
%--------------------------------------------------------------------------
function c   = ctens9(kinematics,properties,dim)
mu1              = 0.595522;
mu2              = 0.050381;
lambda          = properties(3);
J               = kinematics.J;
B               = kinematics.b;
k = 1e5;

for l=1:dim
    for k=1:dim
        for j=1:dim
            for i=1:dim
                sum    =  0;
                eqn = 2*mu2*Bij*Bkl - 1/2*(Bij*Bjl + Bil*Bjk) - ...
                      2/3*(mu1 + 2*mu2*I1)*(Bij*delta(k,l) + Bkl*delta(i,j)) + ...
                      4/3*mu2*(Bij^2*delta(k,l) + Bkl^2*delta(i,j)) + ...
                      2/9*(mu1*I1 + 4*mu2*I2)*delta(i,j)*delta(k,l) + ...
                      1/3*(mu1*I1 + 2*mu2*I2)*(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k));
                c(i,j,k,l) = c(i,j,k,l) + sum;
            end
        end
    end    
end


