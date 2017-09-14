##Sourcing scripts

# source our functions
code.dir <- "c:\temp"
code.files = dir(code.dir, pattern = "[.r]")
for (file in code.files){
  source(file = file.path(code.dir,file))
}