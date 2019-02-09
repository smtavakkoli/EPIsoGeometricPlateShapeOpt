%%
function [Numbering]=numbering(DP1,DP2,NOC1,NOR1,n_e,index)
    Numbering=cell(n_e,1);
    for j=1:n_e

        indj=(mod(j,NOR1)-1+NOR1*logical(~mod(j,NOR1)))+(floor((j-1)/NOR1)*NOC1);

        for i=1:DP2+1
            indi=(i-1)*NOC1;
            for z=1:DP1+1
                ind=indi+indj+z;
                Numbering{j}(i,z)=index(ind);
            end
        end
    end
end