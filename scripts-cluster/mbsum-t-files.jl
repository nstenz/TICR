## julia script to run mbsum to the mb results
## run from the shell like this:
##
## julia path/to/mbsum-t-files.jl mbFolder
##
## The folder named "mbFolder" should have all the tree output files from MrBayes,
## which are named xxx.run1.t, xxx.run2.t etc.
##
## Claudia July 2017, Cecile July 2018

# options to change manually, to fit your analysis:
burnin = 2500 ## 10000*0.25 ## check nex.log file for number of generations etc

length(ARGS) > 0 ||
    error("I need an argument: name of folder containing all the tree sample files")
mbFolder = ARGS[1]

# extract the list of all tree files: xxx.run1.t, xxx.run2.t, etc.
files = String[]
maxNrep = 0
allfiles = filter(x -> endswith(x, ".t"), readdir(mbFolder))
for f in allfiles
    m = match(r"^(.*)\.run([0-9]+)\.t$", f)
    if m!=nothing
        push!(files, m.captures[1])
        nrep = parse(Int, m.captures[2])
        if nrep>maxNrep maxNrep=nrep; end
    end
end
using StatsBase
fc = countmap(files) # file counts: for each locus xxx, count of tree files xxx.run?.t
all([v==maxNrep for v in values(fc)]) ||
    warn("some genes do not have $maxNrep tree files, or their numbers are weird")
info("found $(length(fc)) loci, $maxNrep tree files for each")

for file in keys(fc)
    @show file
    treefiles = join(filter(x -> ismatch(Regex("^$file\.run([0-9]+)\.t\$"), x),
                            allfiles), " ")
    #info("will run: mbsum -n $burnin $treefiles")
    run(`mbsum -n $burnin $treefiles`)
end
