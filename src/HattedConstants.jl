



function hattedConstantsLinear!(k,k̂,u,nMax)
    k̂.=k
    k̂[2] = k̂[2]*sum(u[1:nMax-1])
    k̂[6] = k̂[6]*sum(u[1+nMax:2*nMax-1])
    k̂[10] = k̂[10]*sum(u[1+2*nMax:3*nMax-1])
end

function hattedConstantsNonLinear!(k,k̂,u,nMax)
    k̂.=k
    k̂[2] = k[2]*sum(u[1:nMax-1])
    k̂[3] = k[3]*sum(u[2:nMax].*collect(2:nMax).^(2/3))/sum(u[2:nMax])
    k̂[6] = k[6]*sum(u[1+nMax:2*nMax-1])
    k̂[7] = k[7]*sum(u[nMax+2:2*nMax].*collect(2:nMax).^(2/3))/sum(u[nMax+2:2*nMax])
    k̂[10] = k[10]*sum(u[1+2*nMax:3*nMax-1])
    k̂[11] = k[11]*sum(u[2*nMax+2:3*nMax].*collect(2:nMax).^(2/3))/sum(u[2*nMax+2:3*nMax])
end