function [lambda]=LambdaCalculation(dF,dG,G)
%% Determination of lambda bounds
        b0=dF'*dF;
        lambda=zeros(size(dG,2),1);
    for i=1:size(dG,2)
        a1=dF'*dG(:,i);
        if a1<0
            Bound1=b0/a1;
            a=dG(:,i)'*dG(:,i);
            b=a1+G(i);
            Bound2=b/a;
            lambda(i)=max(Bound1,Bound2)+0.1;
        elseif a1>0
            upperBound=b0/a1;
            a=dG(:,i)'*dG(:,i);
            b=dF'*dG(:,i)+G(i);
            lowerBound=b/a;
            if upperBound>lowerBound
                lambda(i)=(lowerBound+upperBound)/2;
            else
                 lambda(i)=0;
            end
        else
            lambda(i)=0;
        end
    end
end
