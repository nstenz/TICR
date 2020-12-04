#= julia script to paste the mrbayes block in nexus files.
run from the shell like this:

julia path/to/paste-mb-block.jl alignmentFolder

assumptions:
- julia v1 or up (tested on julia v1.5)
- the folder named "alignmentFolder" should have all the alignments.
- alignment files (one file per locus) should have names that
  end with ".nex"
- the current directory should have a MrBayes block, that is,
  a text file named "mb-block.txt" with MrBayes options and commands.
  (else adjust the file name on line 27 below)

The script will create
- a new folder named "alignmentFolder-block"
- new nexus files inside that folder, named: 1.nex, 2.nex, ...
  Such file names make it easier for slurm to use a task-array ID
  as the gene ID (now a number).
- a translate table to map the original names to gene numbers, in
  a new file "translate.txt" (or adjust the file name on line 28)

Claudia 2017-07; Cecile 2018-07, 2020-12
=#

# options to change manually, to fit your file names etc.
mrbayes_block_file = "mb-block.txt"
translate_geneIDs_file = "translate.txt" # this file will be created

length(ARGS) > 0 ||
    error("I need an argument: name of folder containing all the alignments, as nexus files")
nexusFolder = ARGS[1]
nexusFolder = replace(nexusFolder, r"/$" => "") # remove trailing "/" from folder name, if present
@info "will look for nexus files in folder $nexusFolder"
newFolder = string(nexusFolder,"-block")
mkdir(newFolder)
@info "new nexus files will go in $newFolder"

# list of all nexus files
# assumes that the files only have a data block
files = String[]
for f in filter(x -> endswith(x, ".nex"), readdir(nexusFolder))
    push!(files,f)
end
@info "found $(length(files)) nexus files"

# read mrbayes block
block = read(mrbayes_block_file, String)

# translate table file
tab = open(translate_geneIDs_file,"w")
for (i,file) in enumerate(files)
    # @show file
    write(tab, string(i," ",file, "\n"))
    gene = read(joinpath(nexusFolder,file), String)
    newfile = joinpath(newFolder, string(i,".nex"))
    open(newfile,"w") do g
        write(g, gene)
        write(g, block)
    end # closes newfile safely
end
close(tab)

@info "done. check folder $newFolder, and gene IDs in file $translate_geneIDs_file"
