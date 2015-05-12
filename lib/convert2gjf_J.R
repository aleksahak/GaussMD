################################################################################
# Copyright (C) 2012+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the following link:            #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

LAUNCH <- "module add gaussian"
xyz.files <- grep(".xyz",dir(),value=T)

for(file in xyz.files){
  name <- paste(unlist(strsplit(file,".xyz")),".gjf",sep="")
  file.tbl <- read.table(file)
  file.tbl[,2] <- sapply(file.tbl[,2], FUN=function(i){if(i==6){return("C(Iso=13,Spin=1)")};if(i==8){return("O(Iso=17,Spin=5)")};if(i==7){return("N(Iso=15,Spin=1)")};if(i==1){return("H(Iso=1,Spin=1)")};}, USE.NAMES=F, simplify=T)
  file.tbl <- file.tbl[,2:6]
  write.table(file.tbl, file=name, col.names=F, row.names=F, quote=F)

# EXTENSION TO RUN J-COUPLING CALCULATIONS:
  file.lines <- readLines(name)
  file.lines <- c(
  "%nprocshared=1",
  "%mem=1000MB",
  "#P NMR(SpinSpin) SP SCF=Tight RB3LYP/TZVP maxdisk=15GB Prop nosymmetry",
  " ",
  "[From the RB3LYP/6-31G(d,p) MD simulation]",
  " ",
  "0 1",
  file.lines,
  " ")
  
  write(file.lines, file=name)

  rw <- "# you need to submit this with qsub -V scriptname so it picks up"
  rw <- c(rw, paste("# the Gaussian environment variables",sep=""))
  rw <- c(rw, paste(" ",sep=""))
  rw <- c(rw, paste("#PBS -N frm_",unlist(strsplit(name,"_"))[1],sep=""))  #job name
  rw <- c(rw, paste("#PBS -q l1",sep=""))
  rw <- c(rw, paste("#PBS -l walltime=48:00:00",sep=""))
  rw <- c(rw, paste(" ",sep=""))
  rw <- c(rw, paste("# the directory the job was submitted from",sep=""))
  rw <- c(rw, paste("HERE=$PBS_O_WORKDIR",sep=""))
  rw <- c(rw, paste("SCRATCH=${GAUSS_SCRDIR}/${PBS_JOBID}",sep=""))
  rw <- c(rw, paste("mkdir -p $SCRATCH",sep=""))
  rw <- c(rw, paste(" ",sep=""))
  rw <- c(rw, paste("file=",unlist(strsplit(name,".gjf")),sep="")) # filename
  rw <- c(rw, paste("inpfile=${file}.gjf",sep=""))
  rw <- c(rw, paste("outfile=${file}.log",sep=""))
  rw <- c(rw, paste("",sep=""))
  rw <- c(rw, paste("cd $SCRATCH ",sep=""))
  rw <- c(rw, paste("cp $HERE/$inpfile .",sep=""))
  rw <- c(rw, paste(" ",sep=""))
  rw <- c(rw, paste("# Write out some helpful info to the output file",sep=""))
  rw <- c(rw, paste('echo "Starting job $PBS_JOBID"',sep=""))
  rw <- c(rw, paste("echo",sep=""))
  rw <- c(rw, paste('echo "PBS assigned me this node:"',sep=""))
  rw <- c(rw, paste("cat $PBS_NODEFILE",sep=""))
  rw <- c(rw, paste("echo",sep=""))
  rw <- c(rw, paste("",sep=""))
  rw <- c(rw, paste("rm -f $SCRATCH/$outfile",sep=""))
  rw <- c(rw, paste("ln -s $HERE/$outfile $SCRATCH/$outfile ",sep=""))
  rw <- c(rw, paste("",sep=""))
  rw <- c(rw, paste("g03 $inpfile ",sep=""))
  rw <- c(rw, paste("",sep=""))
  rw <- c(rw, paste("rm -f ${SCRATCH}/*",sep=""))
  rw <- c(rw, paste("",sep=""))
  rw <- c(rw, paste("echo",sep=""))
  rw <- c(rw, paste("echo",sep=""))
  rw <- c(rw, paste("qstat -f $PBS_JOBID",sep=""))

  write(rw, file=paste(unlist(strsplit(name,".gjf")),".pbs",sep=""))
  LAUNCH <- c(LAUNCH, paste("qsub -V ",paste(unlist(strsplit(name,".gjf")),".pbs",sep=""),sep=""))
}

write(c(LAUNCH," "), file="launch.scr")
