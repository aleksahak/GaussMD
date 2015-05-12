################################################################################
# Copyright (C) 2012+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the following link:            #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

get.xyz <- function(line.ind=line.ind, md.file=md.file, atom.ctr=DIST.CTR1){
  X <- Y <- Z <- NULL
  for(i in 1:length(atom.ctr)){
    tmp <- linesplit(md.file[line.ind+4+atom.ctr[i]])
    X <- c(X, D2E.numconv(tmp[4]) )
    Y <- c(Y, D2E.numconv(tmp[5]) )
    Z <- c(Z, D2E.numconv(tmp[6]) )
  }
  return(c(mean(X), mean(Y), mean(Z)))
}

