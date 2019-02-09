function [ U ] = CalculateEnergy( P,DIS )

U=0;

for i=1:size(P,2)-1
    
U=U+(P(:,i+1)+P(:,i))'*(DIS(:,i+1)-DIS(:,i))*0.5;

end

end

