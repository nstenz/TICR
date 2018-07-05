## julia script to run mbsum to the mb results
## Claudia July 2017

# list of all nexus files
files = String[]
for f in filter(x -> endswith(x, ".nex"), readdir())
    push!(files,f)
end
println("found $(length(files)) nexus files")
burnin = 2500 ## 10000*0.25 ## from nex.log file

for file in files
    @show file
    f = split(file,".")
    run(`mbsum -n $burnin $(f[1]).nex.run1.t $(f[1]).nex.run2.t $(f[1]).nex.run3.t`)
end

