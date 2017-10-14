using MAT

function extractfrommat(basename,probsizes,algnames)
    iterarrs = Array{Array{Int,2}}(0)
    fevarrs = Array{Array{Int,2}}(0)

    for probname in sort(collect(keys(probsizes)))
        for probsize in probsizes[probname]
            MATfile = matopen(basename*probname*"$(probsize).mat")
            iters = read(MATfile,"iters")
            evals = read(MATfile,"evals")

            push!(iterarrs, hcat(collect(Int.(iters[t]) for t in algnames)...))
            push!(fevarrs, hcat(collect(Int.(evals[t]) for t in algnames)...))
        end
    end
    iters = vcat(iterarrs...)
    fevs = vcat(fevarrs...)
    return iters, fevs
end

function writefile(data, fname,algnames)
    fid = open(fname,"w")
    writedlm(fid, reshape(algnames,1,length(algnames)), ',')
    writedlm(fid, data, ',')
    close(fid)
end

basename = "data/problem"
# TODO: could just loop over all mat-files in data/
#probnames = ["A","B","C","D","E","F","G"]
probsizes = Dict("A" => [100,200], "B" => [100,200], "C" => [100,200],
                 "D" => [500,1000,50000,100000], "E" => [100,200,50000,100000],
                 "F" => [200,500], "G" => [100,200])
algnames = ["ngmreso_sd", "ngmreso_sdls", "ncg", "lbfgs", "ngmres_sd", "ngmres_sdls"]

iters, fevs = extractfrommat(basename,probsizes,algnames)

writefile(iters,"data/iters_AG.csv",algnames)
writefile(fevs,"data/fevs_AG.csv",algnames)


function extractfromcsv(name)
    readdlm(name,',',header=true)[1]
end

#### Create performance profiles

# iters = extractfromcsv("data/iters_AG.csv")
# fevs = extractfromcsv("data/fevs_AG.csv")

using Plots
using BenchmarkProfiles
pgfplots()
#pyplot()
failits = 1500
failidx = find(iters.== failits)

printnames = [ "O-ACCEL-B", "O-ACCEL-A", "NCG", "LBFGS", "N-GMRES-B", "N-GMRES-A"]
iters = convert.(Float64,iters)
fevs = convert.(Float64,fevs)
iters[failidx] .= -1
fevs[failidx] .= -1


(ratios,max_ratio) = performance_ratios(fevs; logscale = false)



# plt1 = performance_profile(iters, printnames, title="Iterations")
# plt2 = performance_profile(fevs, printnames, title="f/g evaluations")

plt = performance_profile(fevs,printnames;title="\$f/g\$ evaluations",
                          logscale=false,
                          sampletol = 3e-2,
                          xscale=:log2,
                          xlims=(1,2^3.5),
                          ylims=(0,1));
algs = [1,5]
pltsd = performance_profile(fevs[:,algs],printnames[algs];title="\$f/g\$ evaluations",
                            logscale=false,
                            sampletol = 2e-2,
                            xscale=:log2, xlims=(1,5),
                            ylims=(0,1));
algs = [2,6]
pltls = performance_profile(fevs[:,algs],printnames[algs],title="\$f/g\$ evaluations";
                            logscale=false,
                            sampletol = 2e-2,
                            xscale=:log2, xlims=(1,6),
                            ylims=(0,1));

savefig(plt,"data/perf_prof_all.tex")
savefig(pltsd,"data/perf_prof_nsd.tex")
savefig(pltls,"data/perf_prof_nls.tex")
