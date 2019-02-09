function [a]=flowvector(sig)
Stress=[sig(1),sig(3),sig(4);sig(3),sig(2),sig(5);sig(4),sig(5),0];
Stress=elhydrpre(Stress);
a1=Stress(1,1);
a2=Stress(2,2);
a3=2*Stress(1,2);
a4=2*Stress(1,3);
a5=2*Stress(2,3);
[~,J2]=vonmises(sig);
a=sqrt(3)*[a1;a2;a3;a4;a5]/(2*sqrt(J2));
end