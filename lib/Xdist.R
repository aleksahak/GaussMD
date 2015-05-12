################################################################################
# Copyright (C) 2012+ Aleksandr B. Sahakyan (aleksahak[at]cantab.net).         #
#                                                                              #
# License: You may redistribute this source code (or its components) and/or    #
#   modify/translate it under the terms of the GNU General Public License      #
#   as published by the Free Software Foundation; version 2 of the License     #
#   (GPL2). You can find the details of GPL2 in the following link:            #
#   https://www.gnu.org/licenses/gpl-2.0.html                                  #
################################################################################

Xdist <- function(Cart1, Cart2) {
 dist <- sqrt( ((Cart1[1] - Cart2[1])^2) +
               ((Cart1[2] - Cart2[2])^2) +
               ((Cart1[3] - Cart2[3])^2)   )
 return(dist)
}