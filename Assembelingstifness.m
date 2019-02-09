function [K]=Assembelingstifness(Ndof,Numbering,KK,DOF)

K=zeros(Ndof);
    %%
    for i=1:size(KK,1)
        numbers=Numbering{i}';
        numbers=numbers(:);
        t=1;
        while ~isempty(numbers)
            R_index=(numbers(1)-1)*DOF+1:numbers(1)*DOF;
            r=(t-1)*DOF+1:t*DOF;
            h=t;
            for j=1:numel(numbers)
                C_index=(numbers(j)-1)*DOF+1:numbers(j)*DOF;
                c=(h-1)*DOF+1:h*DOF;
                if any(r==c)
                    K(R_index,C_index)=K(R_index,C_index)+KK{i}(r,c);
                    else
                    K(R_index,C_index)=K(R_index,C_index)+KK{i}(r,c);
                    K(C_index,R_index)=K(C_index,R_index)+KK{i}(c,r);
                end
                h=h+1;
            end
            t=t+1;
            numbers(1)=[];
        end
    end
end
