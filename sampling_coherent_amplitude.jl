module su2_sampling

include("su2vertex.jl")

using LinearAlgebra
using .su2vertex

export intertwiner_range_result, pre_compute_overlaps, su2_cohAmp_random, thermalization, sampling_int_labels, compute_amp_importance, imp_vertex_amp


function intertwiner_range_result(jays)

    j1,j2,j3,j4 = jays

    min_iota = max(abs(j1-j2),abs(j3-j4))
    max_iota = min(j1+j2,j3+j4)

    return min_iota, max_iota

end


function pre_compute_overlaps(jays,noms) # function that computes the coefficients of coherent intertwiners in intertwiner basis

    j12,j13,j14,j15,j23,j24,j25,j34,j35,j45 = jays

    nl1,nl2,nl3,nl4,nl5 = noms

    #rearrange boundary data from input accordingly
    nl2,nl3,nl4 = nl2[[2,1,4,3]],nl3[[3,4,1,2]],nl4[[4,3,2,1]]

    js = ((j12,j13,j14,j15),(j23,j24,j25,j12),(j34,j35,j13,j23),(j45,j14,j24,j34),(j15,j25,j35,j45))
    vc1 = vector_coherent4j(js[1],nl1)
    vc2 = vector_coherent4j(js[2],nl2)
    vc3 = vector_coherent4j(js[3],nl3)
    vc4 = vector_coherent4j(js[4],nl4)
    vc5 = vector_coherent4j(js[5],nl5)

    range1 = intertwiner_range_result(js[1])
    range2 = intertwiner_range_result(js[2])
    range3 = intertwiner_range_result(js[3])
    range4 = intertwiner_range_result(js[4])
    range5 = intertwiner_range_result(js[5])

    return vc1, vc2, vc3, vc4, vc5, range1, range2, range3, range4, range5

end

#This functions puts out the coherent vectors (vc1,...,vc5) and intertwiner ranges (range1,...,range5) which are the lengths of these vectors


function su2_cohAmp_random(jays,vc1,vc2,vc3,vc4,vc5,range1,range2,range3,range4,range5,Nsample) # function for random sampling coherent vertex amplitude

    j12,j13,j14,j15,j23,j24,j25,j34,j35,j45 = jays

    #vc1 = vector_coherent4j((j12,j13,j14,j15),angs1)
    #vc2 = vector_coherent4j((j23,j24,j25,j12),angs2)
    #vc3 = vector_coherent4j((j34,j35,j13,j23),angs3)
    #vc4 = vector_coherent4j((j45,j14,j24,j34),angs4)
    #vc5 = vector_coherent4j((j15,j25,j35,j45),angs5)

    #range1 = intertwiner_range_result(j12,j13,j14,j15)
    #range2 = intertwiner_range_result(j23,j24,j25,j12)
    #range3 = intertwiner_range_result(j34,j35,j13,j23)
    #range4 = intertwiner_range_result(j45,j14,j24,j34)
    #range5 = intertwiner_range_result(j15,j25,j35,j45)

    result = 0. + 0im

    for samp in 1:Nsample

        i1 = rand(range1[1]:range1[2])
        i2 = rand(range2[1]:range2[2])
        i3 = rand(range3[1]:range3[2])
        i4 = rand(range4[1]:range4[2])
        i5 = rand(range5[1]:range5[2])

        tmp = wig15J(Int(i1),Int(i2),Int(i3),Int(i4),Int(i5),j12,j13,j14,j15,j23,j24,j25,j34,j35,j45) *
        vc1[Int(i1+1-range1[1])]*vc2[Int(i2+1-range2[1])]*vc3[Int(i3+1-range3[1])]*vc4[Int(i4+1-range4[1])]*vc5[Int(i5+1-range5[1])]

        result += tmp

    end

    ranges = length(range1[1]:range1[2]) * length(range2[1]:range2[2]) * length(range3[1]:range3[2]) * length(range4[1]:range4[2]) * length(range5[1]:range5[2])

    result *= ranges / Nsample

    return result

end


function thermalization(vc1,vc2,vc3,vc4,vc5,range1,range2,range3,range4,range5,therm_steps) # function that performs thermalization of the chain, here for all five intertwiner variables

    #vc1 = vector_coherent4j((j12,j13,j14,j15),angs1)
    #vc2 = vector_coherent4j((j23,j24,j25,j12),angs2)
    #vc3 = vector_coherent4j((j34,j35,j13,j23),angs3)
    #vc4 = vector_coherent4j((j45,j14,j24,j34),angs4)
    #vc5 = vector_coherent4j((j15,j25,j35,j45),angs5)

    #range1 = intertwiner_range_result(j12,j13,j14,j15)
    #range2 = intertwiner_range_result(j23,j24,j25,j12)
    #range3 = intertwiner_range_result(j34,j35,j13,j23)
    #range4 = intertwiner_range_result(j45,j14,j24,j34)
    #range5 = intertwiner_range_result(j15,j25,j35,j45)

    count1 = 0

    vcs = [vc1, vc2, vc3, vc4, vc5]

    ranges = [range1 range2 range3 range4 range5]

    iota = [rand(range1[1]:range1[2]) rand(range2[1]:range2[2]) rand(range3[1]:range3[2]) rand(range4[1]:range4[2]) rand(range5[1]:range5[2])]
    lengths = [length(range1[1]:range1[2]) length(range2[1]:range2[2]) length(range3[1]:range3[2]) length(range4[1]:range4[2]) length(range5[1]:range5[2])]
    coh_norm = [sum(abs.(vc1)) sum(abs.(vc2)) sum(abs.(vc3)) sum(abs.(vc4)) sum(abs.(vc5))]

    value_coh4j = [abs(vcs[k][Int(iota[k]+1-ranges[k][1])]) / coh_norm[k] for k in 1:5]

    for j in 1:therm_steps
        rand_i = rand(1:5)
        
        #rand1 = rand(1:max(round(Int,lengths[rand_i]/3),2))
        rand2 = rand(1:2)
        
        iotatest = iota[rand_i] + (-1)^rand2 #* rand1
        if iotatest > ranges[rand_i][2]
            iotatest = ranges[rand_i][2]
        elseif iotatest < ranges[rand_i][1] 
            iotatest = ranges[rand_i][1]
        end

        #println(iotatest)
        
        value_coh4j_new = abs(vcs[rand_i][Int(iotatest+1-ranges[rand_i][1])]) / coh_norm[rand_i]
        r = rand()
        if r <= value_coh4j_new / value_coh4j[rand_i]
            iota[rand_i] = iotatest
            count1 += 1
            value_coh4j[rand_i] = value_coh4j_new
        end
    end

    println("Acceptance rate", count1 / therm_steps)
        
    return iota, vcs, ranges, value_coh4j, coh_norm, lengths
end      


function sampling_int_labels(steps_between,number_samples, iota, value_coh4j, coh_norm, lengths, ranges, vcs) # analogue code for taking samples, parameter to 

    int_samples = zeros(Int,5,number_samples)

    for x in 1:number_samples
    
        for j in 1:steps_between
            rand_i = rand(1:5)
        
            #rand1 = rand(1:max(round(Int,lengths[rand_i]/3),2))
            rand2 = rand(1:2)
        
            iotatest = iota[rand_i] + (-1)^rand2 #* rand1
            if iotatest > ranges[rand_i][2]
                iotatest = ranges[rand_i][2]
            elseif iotatest < ranges[rand_i][1] 
                iotatest = ranges[rand_i][1]
            end

            #println(iotatest)
        
            value_coh4j_new = abs(vcs[rand_i][Int(iotatest+1-ranges[rand_i][1])]) / coh_norm[rand_i]
            r = rand()
            if r <= value_coh4j_new / value_coh4j[rand_i]
                iota[rand_i] = iotatest
                value_coh4j[rand_i] = value_coh4j_new
            end
        end

        int_samples[:,x] = iota

    end

    return int_samples

end


function compute_amp_importance(all_samples,jays,coh_norm,vcs,ranges)  # Given a number of samples, this function computes the coherent vertex amplitude, returns a complex number

    j12,j13,j14,j15,j23,j24,j25,j34,j35,j45 = jays

    result = 0. + 0im
    for samp in 1:length(all_samples[1,:])

        i1,i2,i3,i4,i5 = all_samples[:,samp]

        tmp = wig15J(Int(i1),Int(i2),Int(i3),Int(i4),Int(i5),j12,j13,j14,j15,j23,j24,j25,j34,j35,j45) * prod(coh_norm) *
            vcs[1][Int(i1+1-ranges[1][1])]*vcs[2][Int(i2+1-ranges[2][1])]*vcs[3][Int(i3+1-ranges[3][1])]*vcs[4][Int(i4+1-ranges[4][1])]*vcs[5][Int(i5+1-ranges[5][1])] /
            (abs(vcs[1][Int(i1+1-ranges[1][1])])*abs(vcs[2][Int(i2+1-ranges[2][1])])*abs(vcs[3][Int(i3+1-ranges[3][1])])*abs(vcs[4][Int(i4+1-ranges[4][1])])*abs(vcs[5][Int(i5+1-ranges[5][1])]))

        result += tmp

    end

    result/length(all_samples[1,:])

end


function imp_vertex_amp(jays,vc1,vc2,vc3,vc4,vc5,range1,range2,range3,range4,range5,therm_steps, steps_between, number_samples) # one function to perform thermalization, sampling and calculation of final amplitude

    iota, vcs, ranges, value_coh4j, coh_norm, lengths = thermalization(vc1,vc2,vc3,vc4,vc5,range1,range2,range3,range4,range5,therm_steps)

    all_samples = sampling_int_labels(steps_between,number_samples, iota, value_coh4j, coh_norm, lengths, ranges, vcs)

    result = compute_amp_importance(all_samples,jays,coh_norm,vcs,ranges)

end

end
