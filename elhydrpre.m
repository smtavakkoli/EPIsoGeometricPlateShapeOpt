function [Stress1]=elhydrpre(Stress)

J1=trace(Stress)/3;
Stress1=Stress-J1*eye(3);

end