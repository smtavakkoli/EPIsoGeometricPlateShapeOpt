function [state]=CheckJacobian(detJ)
state=1;
for i=1:length(detJ)
    for j=1:size(detJ{i},1)
        for z=1:size(detJ{i},2)
            if detJ{i}(j,z)<0
                state=0;
                break
        end
    end
end
end