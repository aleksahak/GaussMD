################################################################################
# Copyright (C) 2012+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the following link:            #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

GaussMD <- function(MDOUT="md_start.log", QUERY="traj.query",
                    WIDTH=10, HEIGHT=4, RELATIVE.E=TRUE, CSV=NULL){

# WRITE.QUERY.EXAMPLE = FALSE

###  General script to analyse and plot the ab initio MD output of Gaussian ###
##                                   Supports both separate and joint graphing!
##                                         Supports custom x and y axis ranges!
##                                ?Support graphing from CSV file saved before!
##                       Will support dumping coordinates in different formats!
###  General script to analyse and plot the ab initio MD output of Gaussian ###
source("lib/linesplit.R")
source("lib/get_torsion.R")
source("lib/Xangle.R")
source("lib/Xnormcent.R")
source("lib/D2E_numconv.R")
source("lib/get_xyz.R")
source("lib/Xdist.R")

##-- ** -- ** --
if(is.null(CSV)){
###############################################################################
print("Parsing the query file...", quote=F)

traj.query <- readLines(QUERY)

plot <- grep("PLOT:", traj.query, value=T)
if(length(plot)==0){
  traj.query.plot <- NA
} else {
  write(plot, file="traj.query.plot")
  traj.query.plot <- read.table(file="traj.query.plot", header=T)
  file.remove("traj.query.plot")
}

define <- grep("DEFINE:", traj.query, value=T)
if(length(define)==0){
  traj.query.define <- NA
} else {
  write(define, file="traj.query.define")
  traj.query.define <- read.table(file="traj.query.define", header=T)
  file.remove("traj.query.define")
}

dump <- grep("DUMP:", traj.query, value=T)
if(length(dump)==0){
  traj.query.dump <- NA
} else {
  write(dump, file="traj.query.dump")
  traj.query.dump <- read.table(file="traj.query.dump", header=T)
  file.remove("traj.query.dump")
}
###############################################################################
print("Reading in the Gaussian trajectory...", quote=F)
print("       Make sure you have enough RAM!", quote=F)

md.file <- readLines(MDOUT)
###############################################################################
print("Analysing the Gaussian trajectory...", quote=F)

chk.cart <- grep("Cartesian coordinates:", md.file)
chk.inp  <- grep("  Input orientation:  ", md.file)
chk.distmx <- grep("Distance matrix (angstroms):", md.file, fixed=T)

if(length(chk.inp)>length(chk.cart)){
  chk.inp <- chk.inp[1:length(chk.cart)]
}
Length <- length(chk.cart)
print(paste("There are ",Length," recorded conformations in this trajectory.", sep=""), quote=F)

timestep.line <- grep("Time Step                    =", md.file)
timestep <- as.numeric(linesplit(md.file[timestep.line])[4])

print(paste("The timestep of the simulation is ",timestep," fs.",sep=""), quote=F)
###############################################################################
print("Scanning the trajectory and retrieving time-series data...", quote=F)

data.to.extract <- unique(c(as.character(traj.query.plot[,"Y"]),
                            as.character(traj.query.plot[,"X"]),
                            as.character(traj.query.define[,"GEO.NAME"])))
remove.ind <- c(which(data.to.extract=="Time"), which(data.to.extract=="Etot"),
                which(data.to.extract=="Epot"), which(data.to.extract=="Ekin"),
                which(data.to.extract=="Ekinc"))
if(length(remove.ind)!=0){
 data.to.extract <- data.to.extract[-remove.ind]
}
# Initialising data.to.extract objects
for(i in data.to.extract){
  eval(parse(text=paste(i," <- rep(NA, length=Length)",sep="")))
}

Time  <- seq(from=0, to=(Length*timestep-timestep), by=timestep)
Etot  <- rep(NA, length=Length)
Epot  <- rep(NA, length=Length)
Ekin  <- rep(NA, length=Length)
Ekinc <- rep(NA, length=Length)

DUMP.time   <- DUMP.format <- NULL
######################################
if( (is.na(traj.query.dump[1]))[1] ) {
  print("No molecules are to be extracted from the trajectory.", quote=F)
} else {
  extract.times <- as.character(traj.query.dump[,"TIME"])

  #-- periodic extracts
  periodic.extracts <- grep("per", extract.times, value=T)
  if(length(periodic.extracts)==1){
    periodicity.timestep     <- as.numeric(unlist(strsplit(periodic.extracts,";"))[2])
    periodicity.timestep.adj <- floor(periodicity.timestep/timestep)*timestep
    periodicity.format       <- as.character(traj.query.dump[which(traj.query.dump[,"TIME"]==
                                                                         periodic.extracts ),"FORMAT"])
    print(paste("Geometries will be extracted every ",periodicity.timestep.adj," fs in ",
                                                         periodicity.format," format.",sep=""), quote=F)
    periodicity.times <- seq(from=0, to=(Length*timestep-timestep), by=periodicity.timestep.adj)
    DUMP.time   <- c( DUMP.time, periodicity.times )
    DUMP.format <- c( DUMP.format, rep(periodicity.format, length(periodicity.times)) )

  } else {
    if(length(periodic.extracts)>1){print("ERROR: Only one DUMP line with periodic definition is allowed!", quote=F); break}
  }
  #-- end of periodic extracts

  #-- other extracts
  other.extracts.ind <- which(extract.times!=periodic.extracts)
  if(length(other.extracts.ind)!=0){
    print(paste("Geometries at ",length(other.extracts.ind)," different time-points will be extracted.",sep=""), quote=F)
    DUMP.time   <- c( DUMP.time, floor(as.numeric(as.character(
                                 traj.query.dump[other.extracts.ind,"TIME"]
                                 ))/timestep)*timestep )
    DUMP.format <- c(DUMP.format, as.character(traj.query.dump[other.extracts.ind,"FORMAT"]))
  }
  #-- end of other extracts
  print("  Note, that the extraction step and times can be a bit smaller than the requested ones,", quote=F)
  print("        to make the extraction time-points be multipliers of the MD simulation timestep.", quote=F)
}
######################################


#-- Making the geometry dump folder
if( !is.null(DUMP.time[1]) ){
  system("mkdir geometry_dump")
}

###################
for(i in 1:Length){

  line.ind <- chk.cart[i]
  
  #### Retrieving energy information
  etot <- md.file[(line.ind-20):line.ind]
  etot <- etot[grep(" EKin       =   ", etot)]

  ekin <- as.numeric(unlist(strsplit(linesplit(etot)[3],";")))
  epot <- as.numeric(unlist(strsplit(linesplit(etot)[6],";")))
  etot <- as.numeric(unlist(strsplit(linesplit(etot)[9],";")))

  Etot[i] <- etot
  Ekin[i] <- ekin
  Epot[i] <- epot

  etot <- md.file[(line.ind-20):line.ind]
  etot <- etot[grep(" EKinC      =  ", etot)]
  ekinc <- as.numeric(unlist(strsplit(linesplit(etot)[3],";")))

  Ekinc[i] <- ekinc
  #### 

  ###############################
  df.command.cont <- NULL 
  if(length(data.to.extract)!=0){
  df.command.cont <- paste(", ",paste(paste(data.to.extract,paste("= as.character(",data.to.extract,")",sep="")), collapse=", "), sep="")
  for(geo in data.to.extract){ 

    geo.ind <- NULL
    for(k in 3:6){ 
      geo.ind <- c(geo.ind,
                   as.character( traj.query.define[which(traj.query.define[,"GEO.NAME"]==geo),k]))
    }; geo.ind <- geo.ind[!is.na(geo.ind)]

    geo.ind1       <- as.numeric(unlist(strsplit(geo.ind[1],";")))
    geo.ind1.xyz   <- get.xyz(line.ind=chk.inp[i], md.file=md.file, atom.ctr=geo.ind1)


    if( !is.null(DUMP.time[1]) ){
      #-# GEOMETRY DUMPING
      dump.ind <- which(DUMP.time==Time[i])
      if(length(dump.ind)!=0){
        dump.format <- DUMP.format[dump.ind]
        molecule.extract <- md.file[ (chk.inp[i]+5) :
                        (chk.distmx[which(abs((chk.distmx-chk.inp[i]))==min(abs((chk.distmx-chk.inp[i]))))]-2 )
                            ]
        if(dump.format[1]=="xyz"){write(molecule.extract, file=paste("geometry_dump/",Time[i],"_fs.xyz",sep=""))}
        ### WRITE CASES FOR OTHER FORMATS!!!!
      }
      #-# END OF GEOMETRY DUMPING
    }


    # # # # # # # # # # # #
    if(length(geo.ind)==2){
      # DISTANCE
      geo.ind2     <- as.numeric(unlist(strsplit(geo.ind[2],";")))
      geo.ind2.xyz <- get.xyz(line.ind=chk.inp[i], md.file=md.file, atom.ctr=geo.ind2)      
      measurement  <- Xdist(geo.ind1.xyz, geo.ind2.xyz)
    }
    if(length(geo.ind)==3){
      # ANGLE
      geo.ind2     <- as.numeric(unlist(strsplit(geo.ind[2],";")))
      geo.ind2.xyz <- get.xyz(line.ind=chk.inp[i], md.file=md.file, atom.ctr=geo.ind2) 
      geo.ind3     <- as.numeric(unlist(strsplit(geo.ind[3],";")))
      geo.ind3.xyz <- get.xyz(line.ind=chk.inp[i], md.file=md.file, atom.ctr=geo.ind3)  
      measurement  <- Xangle(Xnormcent(geo.ind2.xyz, geo.ind1.xyz), Xnormcent(geo.ind2.xyz, geo.ind3.xyz))
    }
    if(length(geo.ind)==4){
      # DIHEDRAL ANGLE
      geo.ind2     <- as.numeric(unlist(strsplit(geo.ind[2],";")))
      geo.ind2.xyz <- get.xyz(line.ind=chk.inp[i], md.file=md.file, atom.ctr=geo.ind2) 
      geo.ind3     <- as.numeric(unlist(strsplit(geo.ind[3],";")))
      geo.ind3.xyz <- get.xyz(line.ind=chk.inp[i], md.file=md.file, atom.ctr=geo.ind3) 
      geo.ind4     <- as.numeric(unlist(strsplit(geo.ind[4],";")))
      geo.ind4.xyz <- get.xyz(line.ind=chk.inp[i], md.file=md.file, atom.ctr=geo.ind4)  
      measurement  <- get.torsion(c(geo.ind1.xyz, geo.ind2.xyz, geo.ind3.xyz, geo.ind4.xyz))
    }
    # # # # # # # # # # # #
    eval(parse(text=paste(geo,"[i] <- measurement",sep="")))
  }
  }
  ############################### 

}
###################

df.command <- "results <- data.frame(Time=as.character(Time), Etot=as.character(Etot), Ekin=as.character(Ekin), Ekinc=as.character(Ekinc), Epot=as.character(Epot)" 
df.command <- paste(df.command, df.command.cont, ")",sep="")
eval(parse(text=df.command))
print("Saving the extracted and calculated data into a CSV file.", quote=F)
print("    Note, that the used units in the file are as follows:", quote=F)
print("    energy in a.u (Hartree), angles in degrees, distances in Angstrom, time in fs.", quote=F)
write.csv(results,file=paste(MDOUT,"_result.csv",sep=""), quote=F, row.names=F)
###############################################################################
}
##-- ** -- ** --



###############################################################################
print("Plotting the requested arguments...", quote=F)

results <- read.csv(paste(MDOUT,"_result.csv",sep=""))

plot.groups <- unique(traj.query.plot[,"GROUP"])
for(pg in plot.groups){
  plot.ind <- which(traj.query.plot[,"GROUP"]==pg)

  pdf(width=WIDTH, height=HEIGHT*length(plot.ind),
      file=paste("plot_group_",pg,".pdf", sep=""))
  par(mfrow=c(length(plot.ind),1))
  for(pl in plot.ind){

    X.NAME <- as.character(traj.query.plot[pl,"X"]) 
    Y.NAME <- as.character(traj.query.plot[pl,"Y"])
    X <- results[ ,X.NAME]
#Time     Etot     Ekin    Ekinc      Epot
    if((X.NAME=="Etot"|X.NAME=="Epot"|X.NAME=="Ekin"|X.NAME=="Ekinc")&(RELATIVE.E==TRUE)){
      X <- (X-X[1])*627.503
      print("In the plots, relative energies in kcal/mol are used.", quote=F)
    }
    Y <- results[ ,Y.NAME]
    if((Y.NAME=="Etot"|Y.NAME=="Epot"|Y.NAME=="Ekin"|Y.NAME=="Ekinc")&(RELATIVE.E==TRUE)){
      Y <- (Y-Y[1])*627.503
      print("In the plots, relative energies in kcal/mol are used.", quote=F)
    }

    XLIM=as.character(traj.query.plot[pl,"X.RANGE"])
    if(is.na(XLIM)){
      XLIM <- range(X,na.rm=T)
    } else {
      XLIM <- as.numeric(unlist(strsplit(XLIM,";")))
    }
    YLIM=as.character(traj.query.plot[pl,"Y.RANGE"])
    if(is.na(YLIM)){
      YLIM <- range(Y,na.rm=T)
    } else {
      YLIM <- as.numeric(unlist(strsplit(YLIM,";")))
    }

    plot(x=X, y=Y,
         xlab=as.character(traj.query.plot[pl,"X"]),
         ylab=as.character(traj.query.plot[pl,"Y"]),
         type=as.character(traj.query.plot[pl,"TYPE"]),
         col=as.character(traj.query.plot[pl,"COLOR"]),
         lwd=as.character(traj.query.plot[pl,"LINEWIDTH"]),
         xlim=XLIM, ylim=YLIM)
  }
  dev.off()

}
###############################################################################

}
#pdf(width=7, height=15, file=paste(unlist(strsplit(MDOUT,".log")),".pdf", sep=""))
#  par(mfrow=c(5,1))
#  plot(x=TIME, y=DIST, xlab="Time, fs", ylab="r, A", type="l", col="blue")
#  plot(x=TIME, y=DIST1, xlab="Time, fs", ylab="r1, A", type="l", col="navy")
#  plot(x=TIME, y=EKIN, xlab="Time, fs", ylab="Ekin, a.u.", type="l", col="red")
#  plot(x=TIME, y=(EPOT-min(EPOT))*627.503, xlab="Time, fs", ylab="Epot, kcal/mol", type="l", col="green")
#  plot(x=TIME, y=(ETOT-min(ETOT))*627.503, xlab="Time, fs", ylab="Etot, kcal/mol", type="l", col="black")
#dev.off()

