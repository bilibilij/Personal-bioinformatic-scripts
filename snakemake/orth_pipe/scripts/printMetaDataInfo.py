import subprocess 
import os

def printMeta(config,verion,outfile):
  f = open(outfile, 'w')
  f.write("Synteny pipeline version "+verion+"\n")
  ORTHOFINDER=config['ORTHOFINDER_PROGRAM']
  MACSE=config['MACSE_PROGRAM']
  TREEBEST=config['TREEBEST_PROGRAM']
  IADHORE=config['I-ADHORE_PROGRAM']
  MAFFT=config['MAFFT_PROGRAM']
  bashCommand="\"echo \"Orthofinder version $("+ORTHOFINDER+" --help | grep 'version' | cut -f 3 -d ' ')\" >> "+outfile+"\""
  #print(bashCommand)
  #subprocess.check_output(bashCommand, shell=True)
  bashCommand="echo Orthofinder version $("+ORTHOFINDER+" --help | grep 'version' | cut -f 3 -d ' ')"
  cmd = ['echo Orthofinder']
  with open(outfile, 'w') as out:
    return_code = subprocess.call(cmd, stdout=out)
