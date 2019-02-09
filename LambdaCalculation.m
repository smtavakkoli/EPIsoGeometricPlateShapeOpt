function [lambda]=LambdaCalculation(dF,dG,G)
%% Determination of lambda bounds

    for i=1:size(dG,2)
        a=dF'*dG(:,i);
        b=dF'*dF;
        upperBound=b/a;
        a=dG(:,i)'*dG(:,i);
        b=dF'*dG(:,i)+G;
        lowerBound=b/a;
        lambda(i)=(lowerBound+upperBound)/2;
    end
end