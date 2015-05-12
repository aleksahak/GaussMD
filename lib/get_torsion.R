# ********************************************************************* #
# Borrowed from Bio3D - Barry Grant *C                                  #
# ********************************************************************* #
#########################################################################
## This function calculates  a dihedral angle, given a vector of x, y  ##
##  and z coordiantes of atoms 1, 2, 3 and 4 => c(x1,y1,z1, x2,y2,z2,  ##
##       x3,y3,z3, x4,y4,z4). The returned angle is in degrees.        ## 
##   Getting torsion angle for 1--2--3--4 with respect to 2--3 axis    ##
######################################################################### 
# VERSION 23 Jan, 2010.

get.torsion <- function (xyz, atm.inc = 4) {       # borrowed from bio3d
  if(length(xyz)/3 < 4) {
    return(NULL)
  } else {
    if (!is.vector(xyz) || !is.numeric(xyz)) 
        stop("input 'xyz' should be a numeric vector")
    natm <- length(xyz)/3
    if (natm < 4) 
        stop("Need at least four atoms to define a dihedral")
    if (natm%%1 != 0) 
        stop("There should be three 'xyz' elements per atom")
    m.xyz <- matrix(xyz, nrow = 3)
    atm.inds <- c(1:4)
    out <- NULL
    while (atm.inds[4] <= natm) {
        if (any(is.na(m.xyz[, atm.inds]))) {
            torp <- NA
        }
        else {
            d1 <- m.xyz[, atm.inds[2]] - m.xyz[, atm.inds[1]]
            d2 <- m.xyz[, atm.inds[3]] - m.xyz[, atm.inds[2]]
            d3 <- m.xyz[, atm.inds[4]] - m.xyz[, atm.inds[3]]
            u1 <- (d1[c(2, 3, 1)] * d2[c(3, 1, 2)]) - (d2[c(2, 
                3, 1)] * d1[c(3, 1, 2)])
            u2 <- (d2[c(2, 3, 1)] * d3[c(3, 1, 2)]) - (d3[c(2, 
                3, 1)] * d2[c(3, 1, 2)])
            ctor <- sum(u1 * u2)/sqrt(sum(u1 * u1) * sum(u2 * 
                u2))
            ctor[ctor > 1] <- 1
            ctor[ctor < -1] <- -1
            torp <- matrix(acos(ctor) * (180/pi), ncol = 1)
            if (sum(u1 * ((u2[c(2, 3, 1)] * d2[c(3, 1, 2)]) - 
                (u2[c(3, 1, 2)] * d2[c(2, 3, 1)]))) < 0) 
                torp <- -torp
        }
        out <- c(out, torp)
        atm.inds <- atm.inds + atm.inc
    }
    if (atm.inc == 1 & natm > 4) 
        out <- c(NA, out, NA, NA)
    return(out)
  }
}

#########################################################################
       
