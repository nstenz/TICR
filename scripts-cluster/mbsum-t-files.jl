## julia script to run mbsum on the results of MrBayes
## run from the shell like this:
##
## julia path/to/mbsum-t-files.jl /path/to/mbtrees path/to/mbsumoutput burnin
##
## - The folder "/path/to/mbtrees" should have all the tree output files from MrBayes,
##   which should be named xxx.run1.t, xxx.run2.t etc.
## - The folder "path/to/mbsumoutput" will contain the output of mbsum. Will be created if needed.
## - burnin is optional: 2501 by default.
##
## Assumes: Julia v1 or up; StatsBase package will be added if not done earlier.
##
## Claudia 2017-07, Cecile 2018-07, 2019-07

# read folder names
length(ARGS) >= 2 ||
    error("I need 2 arguments: name of folder containing all the tree sample files, and name of output folder")
mbFolder = ARGS[1]
outfolder = ARGS[2]
isdir(mbFolder)    || error("$mbFolder is not a directory")
!isfile(outfolder) || error("$outfolder is already a file!")
# create the output folder if it doesn't already exist
isdir(outfolder)   || mkdir(outfolder)

# read burnin value, if specified
burnin = 2500+1 ## 10000*0.25 ## check nex.log file for number of generations etc
                # +1 to exclude generation 0
if length(ARGS) >=3
    burnin = parse(Int, ARGS[3])
end

# extract the list of all tree files: xxx.run1.t, xxx.run2.t, etc.
allfiles = filter(x -> endswith(x, ".t"), readdir(mbFolder))
if isempty(allfiles)
    for gene in readdir(mbFolder)
        subfolder = joinpath(mbFolder, gene)
        append!(allfiles, joinpath.(mbFolder, gene, filter(x -> endswith(x, ".t"), readdir(subfolder))))
    end
end
!isempty(allfiles) ||
    error("I could not find the intput tree files in folder $mbFolder nor in its subfolders.")
ntrees = countlines(allfiles[1], eol = '=') - 2 # -1 to exclude final line, -1 to exclude generation 0
@info """Found $ntrees trees in the first .t file.
Will exclude burning=$burnin trees from each .t file.
If that sounds wrong, stop the script with "Control-C" then
re-run with with a second argument to adjust the burnin.
For instance, to exclude 1000 trees from each file, run:
julia $PROGRAM_FILE $mbFolder 1000
"""
maxNrun = 0       # maximum (across genes) number of MrBayes runs per gene
files = String[]  # same list as allfiles, but with ".run?.t" stripped off each file name
for f in allfiles # loop to fill in "files"
    global maxNrun
    m = match(r"^(.*)\.run([0-9]+)\.t$", f)
    if m!=nothing
        push!(files, m.captures[1])
        nruns = parse(Int, m.captures[2])
        if nruns > maxNrun maxNrun=nruns; end
    end
end
# install package StatsBase, if not installed already, and load it
try
    using StatsBase # defines "countmap" used below
catch
    @info "need StatsBase package: will add & use it"
    import Pkg
    Pkg.add("StatsBase")
    using StatsBase
end
# extract unique gene names, and number of runs / gene
fc = countmap(files) # file counts: for each locus xxx, count of tree files xxx.run?.t
all([v==maxNrun for v in values(fc)]) ||
    @warn "some genes do not have $maxNrun tree files, or their numbers are weird"
@info "found $(length(fc)) loci, $maxNrun tree files for each"

# for each gene: run mbsum on all the .run?.t files for that particular gene
re = r"\.run([0-9]+)\.t$"
for gene in keys(fc)
    treefiles = filter(x -> startswith(x, gene) && occursin(re, x), allfiles)
    outfile = joinpath(outfolder, basename(gene) * ".in")
    mbsumcmd = `mbsum -n $burnin -o $outfile $treefiles`
    # @info "will run: $mbsumcmd"
    run(`$mbsumcmd`)
end
