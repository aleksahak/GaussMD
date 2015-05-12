################################################################################
# Copyright (C) 2012+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the following link:            #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

xyz.files <- grep(".xyz",dir(),value=T)

for(file in xyz.files){
  name <- paste(unlist(strsplit(file,".xyz")),".gjf",sep="")
  file.tbl <- read.table(file)
  file.tbl[,2] <- sapply(file.tbl[,2], FUN=function(i){if(i==6){return("C")};if(i==8){return("O")};if(i==7){return("N")};if(i==1){return("H")};}, USE.NAMES=F, simplify=T)
  file.tbl <- file.tbl[,2:6]
  write.table(file.tbl, file=name, col.names=F, row.names=F, quote=F)
}