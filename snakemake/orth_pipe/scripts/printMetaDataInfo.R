outFile <- snakemake@output[[1]]
config <- snakemake@params[["config"]]
cat("Synteny pipeline version",snakemake@params[["verion"]],"\n", file = outFile)
ORTHOFINDER=config[["ORTHOFINDER_PROGRAM"]]
MACSE=config[["MACSE_PROGRAM"]]
TREEBEST=config[["TREEBEST_PROGRAM"]]
IADHORE=config[["I-ADHORE_PROGRAM"]]
MAFFT=config[["MAFFT_PROGRAM"]]
system(paste0("echo Orthofinder version $(",ORTHOFINDER," --help | grep 'version' | cut -f 3 -d ' ') >> ",outFile))
system(paste0("echo MACSE version ",basename(MACSE)," | sed 's#macse_##g' | sed 's#.jar##g' >> ",outFile))
temp_file = tempfile()
system(paste0("echo mafft version $(",MAFFT," --version 2> ",temp_file,"; cat ",temp_file,") >> ",outFile))
system(paste0("echo treebest $(",TREEBEST," 2> ",temp_file," ; grep 'Version' ",temp_file,") >> ",outFile))
system(paste0("echo $(",IADHORE," -v &> ",temp_file,"; grep 'i-ADHoRe' ",temp_file," | sed 's#This is ##g') >> ",outFile))
file.remove(temp_file)

cat("\n=======CONFIG FILE=======\n\n",file = outFile, append = T)
for (n in names(config)){
  if (typeof(config[[n]]) != "list"){
    cat(n,": ",file = outFile,append = T)
    if (length(config[[n]])==1)
      cat(config[[n]],"\n",file = outFile,append = T)  
    else{
      cat("\n",file = outFile,append = T)
      for (ni in config[[n]]){
        cat("-",ni,"\n",file = outFile,append = T)  
      }
    }
  } else {
    cat(n,":\n",file = outFile,append = T)  
    for (nn in names(config[[n]])){
      if (typeof(config[[n]][[nn]]) != "list"){
        cat(" ",config[[n]][[nn]],"\n", file = outFile, append = T)
      }  else {
        cat(" ",nn,":\n",file = outFile,append = T)  
        for (nnn in names(config[[n]][[nn]])){
          cat("  ",nnn,":",config[[n]][[nn]][[nnn]],"\n", file = outFile, append = T)
        }
      }
    } 
  }
}


