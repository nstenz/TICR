## julia script to paste the mrbayes block in nexus files
## you should have in the same directory this script, the mrbayes block
## text file and a folder with the nexus files
## This will name the new nexus files 1.nex, 2.nex,... to make it
## easier for slurm; but it will create a translate table for the original names
## Claudia July 2017

# list of all nexus files
# assumes that the files only have a data block
nexusFolder = "94COS"
files = String[]
for f in filter(x -> endswith(x, "nex"), readdir(nexusFolder))
    push!(files,f)
end
println("found $(length(files)) nexus files")

# read mrbayes block
t = open("mb-block.txt")
block = readstring(t)

# translate table file
tab = open("translate.txt","a+")

i = 1
for file in files
    @show file
    write(tab,string(i," ",file, "\n"))
    f = open(string(nexusFolder,"/",file))
    newfile = string(i,".nex")
    g = open(newfile,"w")
    gene = readstring(f)
    write(g,gene)
    write(g,block)
    close(f)
    close(g)
    i += 1
end
close(tab)

newFolder = string(nexusFolder,"-block")
run(`mkdir $newFolder`)

newfiles = String[]
for f in filter(x -> endswith(x, "nex"), readdir())
    push!(newfiles,f)
end
run(`mv $newfiles $newFolder`)
