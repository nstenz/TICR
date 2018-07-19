## julia script to paste the mrbayes block in nexus files.
## run from the shell like this:
##
## julia path/to/paste-mb-block.jl alignmentFolder
##
## The folder named "alignmentFolder" should have all the alignments.
## The current directory should have a MrBayes block, that is,
## a text file named "mb-block.txt" with MrBayes options and commands.
## The script will create a new folder and new nexus files named:
## 1.nex, 2.nex, ... These new file names will make it easier for slurm
## to use the slurm task-array ID as the gene ID (now a number).
## A translate table is create to map the original names to gene numbers
## Claudia July 2017 - Cecile July 2018

# options to change manually, to fit your file names etc.
mrbayes_block_file = "mb-block.txt"
translate_geneIDs_file = "translate.txt" # this file will be created

length(ARGS) > 0 ||
    error("I need an argument: name of folder containing all the alignments, as nexus files")
nexusFolder = ARGS[1]
nexusFolder = replace(nexusFolder, r"/$", "") # remove trailing "/" from folder name, if present
info("will look for nexus files in folder $nexusFolder")
newFolder = string(nexusFolder,"-block")
mkdir(newFolder)
info("new nexus files will go in $newFolder")

# list of all nexus files
# assumes that the files only have a data block
files = String[]
for f in filter(x -> endswith(x, "nex"), readdir(nexusFolder))
    push!(files,f)
end
info("found $(length(files)) nexus files")

# read mrbayes block
t = open(mrbayes_block_file)
block = readstring(t)
close(t)

# translate table file
tab = open(translate_geneIDs_file,"w")
i = 1
for file in files
    # @show file
    write(tab, string(i," ",file, "\n"))
    open(joinpath(nexusFolder,file)) do f
        gene = readstring(f)
        newfile = joinpath(newFolder, string(i,".nex"))
        open(newfile,"w") do g
          write(g, gene)
          write(g, block)
        end # closes newfile safely
    end # closes f safely
    i += 1
end
close(tab)

info("done. check folder $newFolder, and gene IDs in file $translate_geneIDs_file")
