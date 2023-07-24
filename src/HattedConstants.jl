



function hattedConstants!(k,k̂,u,nMax)
    k̂.=k
    k̂[2] = k̂[2]*sum(u[1:nMax-1])
    k̂[6] = k̂[6]*sum(u[1+nMax:2*nMax-1])
    k̂[10] = k̂[10]*sum(u[1+2*nMax:3*nMax-1])
end