%--------------------------------------------------------------------------
% Evaluates the Cauchy stress tensor for material type 8.
%--------------------------------------------------------------------------
function Cauchy = stress9(kinematics,properties,cons)
mu1              = 0.595522;
mu2              = 0.050381;
lambda          = properties(3);
J               = kinematics.J;
b               = kinematics.b;
k = 1e5;

% Additions for Mooney Rivilin

C               = kinematics.F'*kinematics.F;
C_bar           = (J^(-2/3)*C);

Siso            = J^(-2/3) * (((-1/3) * (mu1 * cons.I + 2 * mu2 * cons.I) * (C_bar^-1)) + ((mu1 + mu2 * cons.I) * cons.I) - (mu2 * C_bar));
Svol            = (J^(1/3)) * (J - 1) * k * (C_bar^-1)
S               = Siso + Svol
Cauchy          = kinematics.F*S*kinematics.F';

end


