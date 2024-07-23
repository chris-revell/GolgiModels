#
#  CompartmentGeometry.jl
#  GolgiModels
#
#  Created by Christopher Revell on 13/12/2023.

function ellipsoidAxes(cellType,maturation,V)
    if maturation == "Trans"
        r=0.424
    elseif maturation == "Medial"
        r=0.406
    elseif maturation == "Cis"
        r=0.419
    end
    a = cbrt(3.0*V/(4.0*π*r))
    b = r*a    
    return (a,b)
end

function ellipsoidEquator(a,b)
    h = ((a-b)^2)/((a+b)^2)
    c = π*(a+b)*(1+3.0*h/(10+sqrt(4-3.0*h)))
    return c
end
