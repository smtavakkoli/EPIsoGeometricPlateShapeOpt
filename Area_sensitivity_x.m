function [A_px]= Area_sensitivity_x(wanted,Re,Rn,Xe,Xn,Ye,Yn,...
    detJ2,Numbering,w_guss_x,w_guss_y)
%%
A_px=zeros(length(wanted),1);

  for ii=1:length(wanted)
      C=[];
for i=1:length(Numbering)
    ind=find(Numbering{i}(:)==wanted(ii));
    if ~isempty(ind)
        C=[C;i];
    end
end
  
C(isempty(C))=[];
C(C==0)=[];
    %% Calculate K_px{i} and K_py{i}
    Re1=Re(C);Rn1=Rn(C);
    
    Ye1=Ye(C);Yn1=Yn(C);
    Numbering1=Numbering(C);

    detJ2_1=detJ2(C);
    wgx=w_guss_x(C);
    wgy=w_guss_y(C);
    %%...........................
    for i1=1:length(C)
        
        [r,c]=find(Numbering1{i1}==wanted(ii));

   %%................................     
        for countx=1:length(wgx{i1})
            for county=1:length(wgy{i1})
    %%..................................

                Xe_px=Re1{i1}{countx,county}(r,c);
                Xn_px=Rn1{i1}{countx,county}(r,c);
           
    %%..............................................................................
yn1=Yn1{i1}(countx,county);
ye1=Ye1{i1}(countx,county);

                detJ1_px=Xe_px*yn1-Xn_px*ye1;
                
   %%.................................................................................

 A_px(ii,1)=A_px(ii,1)+detJ1_px*wgx{i1}(countx)*wgy{i1}(county)*detJ2_1(i1);

  %%.................................................................................              
            end

        end

    end
  end
end