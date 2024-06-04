using LinearAlgebra
using SparseArrays

# A is a matrix of nVerts x nEdges
# A[1:nEdgesi, :] corresponds to i directed edges
# A[1:nEdgesj, :] corresponds to i directed edges
# A[1:nEdgesk, :] corresponds to i directed edges
# Vertices and edges in each direction ordered in 1D vector as per julia system data layout

ni=13
nj=13
nk=1

nVerts = ni*nj*nk
nEdgesi = (ni-1)*nj*nk
nEdgesj = ni*(nj-1)*nk
nEdgesk = ni*nj*(nk-1)
nEdges = nEdgesi+nEdgesj+nEdgesk

A = spzeros(Int64, nEdges, nVerts)
vertArray = reshape(collect(1:nVerts), (ni, nj, nk))
iEdgeArray = reshape(collect(1:nEdgesi), (ni-1, nj, nk))
jEdgeArray = reshape(collect(1:nEdgesj), (ni, nj-1, nk))
kEdgeArray = reshape(collect(1:nEdgesk), (ni, nj, nk-1))

# Link i-directed edges to corresponding vertices
for kk=1:nk
    for jj=1:nj
        for ii=1:(ni-1)
            edgeIndex = 1 + (ii-1)*stride(iEdgeArray,1) + (jj-1)*stride(iEdgeArray,2) + (kk-1)*stride(iEdgeArray,3)
            v1 = 1 + (ii-1)*stride(vertArray,1) + (jj-1)*stride(vertArray,2) + (kk-1)*stride(vertArray,3)
            v2 = 1 + (ii)*stride(vertArray,1) + (jj-1)*stride(vertArray,2) + (kk-1)*stride(vertArray,3)
            A[edgeIndex,v1] = -1
            A[edgeIndex,v2] = 1
        end
    end
end  
# Link j-directed edges to corresponding vertices
for kk=1:nk
    for jj=1:(nj-1)
        for ii=1:ni
            edgeIndex = nEdgesi + 1 + (ii-1)*stride(jEdgeArray,1) + (jj-1)*stride(jEdgeArray,2) + (kk-1)*stride(jEdgeArray,3)
            v1 = 1 + (ii-1)*stride(vertArray,1) + (jj-1)*stride(vertArray,2) + (kk-1)*stride(vertArray,3)
            v2 = 1 + (ii-1)*stride(vertArray,1) + (jj)*stride(vertArray,2) + (kk-1)*stride(vertArray,3)
            A[edgeIndex,v1] = -1
            A[edgeIndex,v2] = 1
        end
    end
end  
# Link k-directed edges to corresponding vertices
for kk=1:(nk-1)
    for jj=1:nj
        for ii=1:ni
            edgeIndex = nEdgesi + nEdgesj + 1 + (ii-1)*stride(kEdgeArray,1) + (jj-1)*stride(kEdgeArray,2) + (kk-1)*stride(kEdgeArray,3)
            v1 = 1 + (ii-1)*stride(vertArray,1) + (jj-1)*stride(vertArray,2) + (kk-1)*stride(vertArray,3)
            v2 = 1 + (ii-1)*stride(vertArray,1) + (jj-1)*stride(vertArray,2) + (kk)*stride(vertArray,3)
            A[edgeIndex,v1] = -1
            A[edgeIndex,v2] = 1
        end
    end
end  
