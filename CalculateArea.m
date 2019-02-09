function [PlateArea]=CalculateArea(detJ1,detJ2,guss_x,w_guss_x,guss_y,w_guss_y)

PlateArea=0;

for i=1:length(detJ1)

for countx=1:length(guss_x{i})
    
    for county=1:length(guss_y{i})
        
    PlateArea=PlateArea+detJ1{i}(countx,county)*detJ2(i)*w_guss_x{i}(countx)*w_guss_y{i}(county);

    end
end
end