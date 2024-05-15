#
#  MakeIncidenceMatrix.jl
#  GolgiModels
#
#  Created by Christopher Revell on 02/05/2024.

module MakeIncidenceMatrix

using LinearAlgebra
using SparseArrays

# State matrix uMat has dimensions Nplus x Nplus. 
# Column major order, so when flattened to a state vector u, u[1:Nplus] == uMat[:,1], u[1+Nplus:2*Nplus]==uMat[2,:] etc
# State vector has size Nplus*Nplus
# Let us say that in the state matrix, the first dimension is ν and the second is x.
# So in the state vector, the first Nplus components correspond to 0<=ν<=1.0 for x=0.0; second Nplus components correspond to 0<=ν<=1.0 for x=dx, and so on 

# Function to convert index (i,j) to access uMat[i,j] where uMat is of size N1xN2 to index k, 
# which is the corresponding index of the same data point when uMat is flattened to a vector 
ij_To_k(i, j, N1) = (j-1)*N1+i 
# Function to convert index k to access u[k] where u is a vector of size N1*N2 to index (i,j) to access the same data 
# in uMat, where uMat is u reshaped to be a 2D matrix of size N1xN2 
function k_To_ij(k, N1)
    j = floor(Int64, (k-1)/N1)
    i = k-j*N1
    return (i,j)
end

# https://mbernste.github.io/posts/laplacian_matrix/
# Make incidence matrix mapping grid points to directed edges connecting them 
# Edges oriented from low index towards high index.
# Edges labelled such that all i-directed edges form the first half of the edge list, in column major order. 
function makeIncidenceMatrix2D(Nplus)

    #      jj->
    #    *---*---*
    # ii \   \   \
    # \  *---*---*
    # v  \   \   \
    #    *---*---*
 
    nVerts = Nplus*Nplus
    edgesPerDirection = (Nplus-1)*Nplus
    nEdges = 2*edgesPerDirection
    dimEdges1 = Nplus-1
    dimEdges2 = Nplus
    A = spzeros(Int64, nEdges, nVerts)
    # 1st dimension i directed Edges
    for ii=1:dimEdges1
        for jj=1:dimEdges2
            j  = ij_To_k(ii, jj, dimEdges1)
            k1 = ij_To_k(ii, jj, Nplus)
            k2 = ij_To_k(ii + 1, jj, Nplus)
            A[j, k1] = -1
            A[j, k2] = 1
        end
    end 
    # 2nd dimension j directed Edges
    for ii=1:Nplus
        for jj=1:dimEdges1
            j  = ij_To_k(ii, jj, Nplus) + edgesPerDirection
            k1 = ij_To_k(ii, jj, Nplus)
            k2 = ij_To_k(ii, jj + 1, Nplus)
            A[j, k1] = -1
            A[j, k2] = 1
        end
    end 
    return A
end


# A is a matrix of nVerts x nEdges
# A[1:nEdgesi, :] corresponds to i directed edges
# A[1:nEdgesj, :] corresponds to i directed edges
# A[1:nEdgesk, :] corresponds to i directed edges
# Vertices and edges in each direction ordered in 1D vector as per julia system data layout
function makeIncidenceMatrix3D(ni, nj, nk)

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

    return A
end

export ij_To_k, k_To_ij, makeIncidenceMatrix2D, makeIncidenceMatrix3D

end