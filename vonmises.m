function [sigeffect,J2]=vonmises(sig)

Stress=[sig(1),sig(3),sig(4);sig(3),sig(2),sig(5);sig(4),sig(5),0];
[Stress]=elhydrpre(Stress);
J2=0.5*(Stress(1,1)^2+Stress(2,2)^2+Stress(3,3)^2+...
Stress(1,2)^2+Stress(2,1)^2+Stress(1,3)^2+...
Stress(3,1)^2+Stress(3,2)^2+Stress(2,3)^2);
sigeffect=sqrt(3*J2);

end