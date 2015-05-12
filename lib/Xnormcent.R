################################################################################
# Copyright (C) 2012+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the following link:            #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

Xnormcent <- function(Vec1B, Vec1E) {
 L.x <- Vec1E[1] - Vec1B[1]
 L.y <- Vec1E[2] - Vec1B[2]
 L.z <- Vec1E[3] - Vec1B[3]
 L <- Xdist(Vec1B, Vec1E)
 normcent.xyz <- c(L.x/L, L.y/L, L.z/L)
 return(normcent.xyz)
}

################################################################################