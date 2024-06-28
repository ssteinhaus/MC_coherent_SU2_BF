module su2vertex

using LinearAlgebra, WignerSymbols, HalfIntegers, Combinatorics, TensorOperations, Memoization

export intw_range, dim_j, wignermatrix, coherent4j, vector_coherent4j, wig15J, tensor_15j, cohn_vertex

#========= Utilities  ===============# 
# check integerness and (j,m) as angular momentum
isAng(j, m) = (abs(m) <= j && ishalfinteger(j) && isinteger(j-m) && isinteger(j+m))

# check integerness and (j,m) as angular momentum
isPair(j1, j2) = (ishalfinteger(j1) && isinteger(j1-j2) && isinteger(j1+j2))

# intertwiner range 
intw_range(j1,j2,j3,j4) =  max(abs(j1-j2),abs(j3-j4)):min(j1+j2,j3+j4)


# intertwiner range 
intw_range(j1,j2,j3,j4,j5,j6) =  max(abs(j1-j2),abs(j3-j4),abs(j5-j6)):min(j1+j2,j3+j4,j5+j6)

# spin dimension 
dim_j(j) = twice(j)+1.0

#================= Coherent 4j amplitude ================#
#get s2 angles (ϕ,θ) from a 3d (unit) normal vector
function s2angles(n)# input is a 3D unit normal vector
    n = normalize(n)
    return [atan(n[2],n[1]),acos(n[3])] # give conditions for n[1]==0
end

# wigner matrix as function of S^2 angles (ϕ,θ).. |+> sector  
function wignermatrix_ang(j,m,phi,theta)::ComplexF64 #wigner matrix as function of angles on S^2.. gives g|+>
    j <= 33 ? two_j = Int(2*j) : two_j = big(Int(2*j))
    return sqrt(binomial(two_j,Int(j+m)))* (cos(theta/2))^(j+m)*(sin(theta/2)*exp(-im*phi))^(j-m) 
end

# wigner matrix as function of 3d (unit) normal vector.. |+> sector 
function wignermatrix(j,m,nv)::ComplexF64
   return wignermatrix_ang(j,m,atan(nv[2],nv[1]),acos(nv[3]))
end


#4j intertwiner at channel j12 - for all incoming edges
@memoize function wig3j(j1,j2,j3,m1,m2,m3=-m1-m2)::Float64
    wigner3j(j1,j2,j3,m1,m2,m3)
end

function wigner4j(i,j1,j2,j3,j4,m1,m2,m3,m4=-m1-m2-m3)::Float64
    isAng(i,m1+m2) ? (-1)^(i+m1+m2)*wig3j(j1,j2,i,m1,m2,-m1-m2)*wig3j(i,j3,j4,m1+m2,m3,m4) : 0.0 
end

function coherent4j(iota,jay,nvecs)::ComplexF64
    r::ComplexF64 = 0.0
    j1,j2,j3,j4 = jay
    n1,n2,n3,n4 = nvecs
    for m1 in -j1:j1
        wm1 = wignermatrix(j1,m1,n1)
        if wm1 != 0.0
            for m2 in -j2:j2
                wm2 = wignermatrix(j2,m2,n2)
                if wm2 != 0.0 
                    @simd for m3 in max(-j3, -m1-m2-j4):min(j3,j4-m1-m2)
                        #m4 = -(m1+m2+m3)
                        wm34 = wignermatrix(j3,m3,n3)*wignermatrix(j4,-(m1+m2+m3),n4)
                        w4j = wigner4j(iota,j1,j2,j3,j4,m1,m2,m3,)
                        if wm34 != 0.0 && w4j != 0.0
                            r += wm1*wm2*wm34*w4j
                        end
                    end
                end
            end
        end
    end
    return r
end

#function vector_coherent4j(js,nvecs)::Vector{ComplexF64}
#    [dim_j(i)*coherent4j(i,js,nvecs) for i in intw_range(js[1],js[2],js[3],js[4])]
#end

function vector_coherent4j(js,nvecs)::Vector{ComplexF64}
    int_labels = intw_range(js[1],js[2],js[3],js[4])

    vec_coh = complex(zeros(length(int_labels)))

    Threads.@threads for i in int_labels

        vec_coh[Int(i - int_labels[1]) + 1] += dim_j(i)*coherent4j(i,js,nvecs)

    end

    return vec_coh
    
end


#----------------- Wigner 6J and 15J symbols  -------------------------------
# 6J symbol 
@memoize function wig6j(j1,j2,j3,j4,j5,j6)::Float64
#function wig6j(j1,j2,j3,j4,j5,j6)::Float64
    wigner6j(j1,j2,j3,j4,j5,j6)
end

# 15J symbol 
function wig15J(i1,i2,i3,i4,i5,j12,j13,j14,j15,j23,j24,j25,j34,j35,j45)::Float64
    sgn = i1+i2+i3+i4+i5+j12+j13+j14+j15+j23+j24+j25+j34+j35+j45
    sol = 0.0
    if isinteger(sgn)
        ss_min = max(max(max(max(abs(i1-j25),abs(i2-j13)),abs(i3-j24)),abs(i4-j35)),abs(i5-j14))
        ss_max = min(min(min(min(i1+j25,i2+j13),i3+j24),i4+j35),i5+j14)
        @simd for x in ss_min:ss_max 
            @inbounds sol += dim_j(x)*wig6j(i1,j25,x,i2,j13,j12)*wig6j(i2,j13,x,i3,j24,j23)*
            wig6j(i3,j24,x,i4,j35,j34)*wig6j(i4,j35,x,i5,j14,j45)*wig6j(i5,j14,x,i1,j25,j15)
        end
        return (-1)^sgn *sol 
    else
        sol
    end
end

# 15J symbol as tensor -- note CONVENTION can be read from the order of spins in intertwiners
function tensor_15j(j12,j13,j14,j15,j23,j24,j25,j34,j35,j45)
    #js = ((j12,j13,j14,j15),(j23,j24,j25,j12),(j34,j35,j13,j23),(j45,j14,j24,j34),(j15,j25,j35,j45))
    i1,i2,i3,i4,i5 = intw_range(j12,j13,j14,j15),intw_range(j23,j24,j25,j12),intw_range(j34,j35,j13,j23),
    intw_range(j45,j14,j24,j34),intw_range(j15,j25,j35,j45)
    sol = Array{Float64}(undef,length.((i1,i2,i3,i4,i5)))
    for (k1,v1) in enumerate(i1) 
        for (k2,v2) in enumerate(i2) 
            for (k3,v3) in enumerate(i3) 
                for (k4,v4) in enumerate(i4) 
                    for (k5,v5) in enumerate(i5)
                        @inbounds sol[k1,k2,k3,k4,k5] = 
                        wig15J(v1,v2,v3,v4,v5,j12,j13,j14,j15,j23,j24,j25,j34,j35,j45)
                    end
                end
            end
        end
    end
    sol
end

function cohn_vertex(jays,noms) # note positions of noms and jays
    # lst is a list of which edges to which noms connect consistently
    #@assert length(noms) == length(lens)
    j12,j13,j14,j15,j23,j24,j25,j34,j35,j45 = jays
    nl1,nl2,nl3,nl4,nl5 = noms
    #rearrange boundary data from input accordingly
    nl2,nl3,nl4 = nl2[[2,1,4,3]],nl3[[3,4,1,2]],nl4[[4,3,2,1]]
    ct = tensor_15j(j12,j13,j14,j15,j23,j24,j25,j34,j35,j45)
    #tnsr15 = tensor_15j(j12,j13,j14,j15,j23,j24,j25,j34,j35,j45)
    js = ((j12,j13,j14,j15),(j23,j24,j25,j12),(j34,j35,j13,j23),(j45,j14,j24,j34),(j15,j25,j35,j45))
    vc1 = vector_coherent4j(js[1],nl1)
    vc2 = vector_coherent4j(js[2],nl2)
    vc3 = vector_coherent4j(js[3],nl3)
    vc4 = vector_coherent4j(js[4],nl4)
    vc5 = vector_coherent4j(js[5],nl5)
    @tensor ct[a,b,c,d,e]*vc1[a]*vc2[b]*vc3[c]*vc4[d]*vc5[e]
end



end