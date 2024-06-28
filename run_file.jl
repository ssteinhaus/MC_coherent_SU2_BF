include("su2vertex.jl")

include("sampling_coherent_amplitude.jl")

using Plots, BenchmarkTools, JLD2, Memoization
using .su2vertex
using .su2_sampling

#Equilateral 4-simplex

#nn= [[0.0, 0.0, 1.0], [0.0, 0.9428090415820634, -0.3333333333333334], 
#    [0.816496580927726, -0.4714045207910318, -0.3333333333333332], 
#    [-0.816496580927726, -0.47140452079103157, -0.3333333333333332]]

#Full calculation

#xax= 20.5:0.5:30
#amp1= []
#@time for i in xax
#    push!(amp1,cohn_vertex(i*ones(10),[nn,nn,nn,nn,nn]) )
#    Memoization.empty_all_caches!();
#end

#save_object("results/equi_full_results_30", amp1);

#MC calculation

#xax_imp= 39.5:0.5:40
#amp_imp_equi= complex(zeros(length(xax_imp),30))
#for i in xax_imp
#    @time vc1, vc2, vc3, vc4, vc5, range1, range2, range3, range4, range5 = pre_compute_overlaps(i*ones(10),[nn,nn,nn,nn,nn])
#
#    for l in 1:30
#
#        @time amp_imp_equi[Int(2*i + 1)-79, l] =  imp_vertex_amp(i*ones(10),vc1,vc2,vc3,vc4,vc5,range1,range2,range3,range4,range5,10^4, 1000, 10^5) 
#
#    end
#    Memoization.empty_all_caches!();
#end

#save_object("results/test_equi_mc_results_1e5_100", amp_imp_equi);


#Isosceles 4-simplex
#boundary data for tetrahedra
# equilateral [1,1,1,1]
n1 = [[0.0, 0.0, 1.0],
 [0.0, 0.9428090415820634, -0.33333333333333326],
 [0.816496580927726, -0.4714045207910317, -0.3333333333333335],
 [-0.816496580927726, -0.4714045207910317, -0.33333333333333337]]

# isosceles [1,1,1,2]
n2 = [[0.0, 0.0, 1.0],
 [0.0, 0.9860132971832694, 0.16666666666666682],
 [0.9759000729485332, 0.14085904245475284, 0.16666666666666682],
 [-0.4879500364742665, -0.563436169819011, -0.6666666666666667]];


 ns1,ns2,ns3,ns4,ns5 = [n2,n2,n2,n2,n1];


#amp = []
#bdyars = [1,1,1,2,1,1,2,1,2,2]
#@time for i in 31:35
#    push!(amp,cohn_vertex(i*bdyars,[ns1,ns2,ns3,ns4,ns5]) )
#end


#save_object("results/isosceles_full_results_35", amp);

bdyars = [1,1,1,2,1,1,2,1,2,2]
xax_imp= 116:120
amp_imp_iso= complex(zeros(length(xax_imp),30))
println(xax_imp);
@time for i in xax_imp
    vc1, vc2, vc3, vc4, vc5, range1, range2, range3, range4, range5 = pre_compute_overlaps(i*bdyars,[ns1,ns2,ns3,ns4,ns5])

    Threads.@threads for l in 1:30

        amp_imp_iso[i+1 - 116, l] =  imp_vertex_amp(i*bdyars,vc1,vc2,vc3,vc4,vc5,range1,range2,range3,range4,range5,10^4, 1000, 10^5)

    end
    Memoization.empty_all_caches!();
end

save_object("results/isosceles_mc_results_1e5_120", amp_imp_iso);


#Non-regular 4-simplex
#boundary data for tetrahedra

#spins [2,1,1,1]
#nv1 = [[0.0, 0.0, 1.0], [0.0, 0.830990198869012, -0.5562870566386028], 
#    [0.44286684462550546, -0.7031455528891631, -0.5562870566386043], 
#    [-0.4428668446255052, -0.12784464597984913, -0.8874258867227933]];
#spins [2,2,2,1]
#nv2 = [[0.0, 0.0, 1.0], [0.0, 0.8887803753208977, -0.4583333333333331], 
#    [0.4736654667156709, -0.7520449329638362, -0.45833333333333387], 
#    [-0.9473309334313413, -0.273470884714125, -0.16666666666666644]];
#spins [1,2,1,1]
#nv3 =  [[0.0, 0.0, 1.0], [0.0, 0.8309901988690119, -0.556287056638603], 
#    [0.4428668446255053, -0.8534198618296067, -0.274851773445586], 
#    [-0.442866844625505, -0.8085605359084173, 0.3874258867227926]];
#spins [1,2,1,1]
#nv4 = [[0.0, 0.0, 1.0], [0.0, 0.8309901988690113, -0.5562870566386039], 
#    [0.4428668446255061, -0.8534198618296062, -0.27485177344558676], 
#    [-0.44286684462550624, -0.808560535908416, 0.3874258867227939]];
#spins [1,1,1,1]
#nv5 = [[0.0, 0.0, 1.0], [0.0, 0.9428090415820635, -0.3333333333333334], 
#    [0.816496580927726, -0.47140452079103295, -0.3333333333333318], 
#    [-0.8164965809277261, -0.4714045207910306, -0.3333333333333348]];



#amp = [];
#bdyars = [2,1,1,1,2,2,1,1,1,1];
#@time for i in 21:30
#    push!(amp,cohn_vertex(i*bdyars,[nv1,nv2,nv3,nv4,nv5]) )
#end

#save_object("results/non_reg_full_results_30", amp);



#bdyars = [2,1,1,1,2,2,1,1,1,1];
#xax_imp= 121:130
#amp_imp= complex(zeros(length(xax_imp),30))
#println(xax_imp);
#@time for i in xax_imp
#    vc1, vc2, vc3, vc4, vc5, range1, range2, range3, range4, range5 = pre_compute_overlaps(i*bdyars,[nv1,nv2,nv3,nv4,nv5])

#    Threads.@threads for l in 1:30

#        amp_imp[i+1 - 121, l] =  imp_vertex_amp(i*bdyars,vc1,vc2,vc3,vc4,vc5,range1,range2,range3,range4,range5,10^4, 1000, 10^5)

#    end
#    Memoization.empty_all_caches!();
#end

#save_object("results/non_reg_mc_results_1e5_130", amp_imp);