# Computation of "Effort-To-Compress" (ETC) measure using Non-Sequential Recursive Pair
# Substitution (NSRPS)
#
# ETC <- function (InputTS, ...)
# 
# InputTS - input time series of real numbers or a symbolic sequence of integers
# NumBins - Number of Bins to partition the input time series to create a symbolic sequence,
#           if this is set to 0 then input is already a symbolic sequence
#
# OUTPUTS:
#   [[1]] N - ETC measure = number of iterations for entropy -> 0 in the NSRPS algorithm
#   [[2]] EntLen_array - Array of entropy x length values for the NSRPS algorithm    
#   
#   Example1:
#     > ETC(c(1, 2, 1, 2, 2, 2, 2), 0)
#     gives:
#     [[1]] [1] 5
#     [[2]] [1] 6.041844 4.854753 6.000000 2.754888 2.000000 0.000000
# 
#   Example2:
#     > ETC(seq(0, .99, .1), 2)
#     gives:
#     [[1]] [1]
#     [[2]] [1] 10.000000 10.390360 11.509775  9.609640  8.000000  4.754888  2.000000  0.000000
#       
#     by: Nithin Nagaraj, Feb. 2013
#     Amrita Vishwa Vidyapeetham, Amritapuri, Kollam, Kerala, INDIA
#       
#     Ref:
#     1. Nithin Nagaraj, Karthi Balasubramanian, Sutirth Dey, A New Complexity Measure For 
#     Time Series Analysis and Classification, Eur. Phys. J. Special Topics 222, 847-860 (2013)

#     2. W. Ebeling, M.A. Jiménez-Montano, On Grammars, Complexity, and 
#     Information Measures of Biological Macromolecules, Mathematical Biosciences 52, (1980) 53-71.
#
# The latest version of this program can be downloaded from the following websites:
# 
# URL1: https://sites.google.com/site/nithinnagaraj2/journal/etc
# URL2: https://sites.google.com/a/acads.iiserpune.ac.in/sdlab/publications?pli=1
# 
# This program may be used for research purposes only. We provide no guarantees on
# the validity of the output, or the interpretation of results.
# Please do not remove this notice. 
#
# v1.0 (Feb 15, 2013)
# v1.1 (May 31, 2013)
#
# Transcription from Matlab: Daniel Koska, Chemnitz University of Technology, Feb. 2021.

# ------------------------------------------------------------------------------
# Wrapper function
# ------------------------------------------------------------------------------
ETC <- function (InputTS, ...) {
  
  # What am I doing here?
  args <- list(...)
  
  if (length(args) == 0) {
    NumBins <- 4
  } else if (is.atomic(args[[1]])
             && length(args[[1]]) == 1L
             && round(args[[1]]) == args[[1]]
             && is.numeric(args[[1]])
             && args[[1]] >= 0) {
    NumBins <- args[[1]]
  } else {
    stop('NumBins should an integer greater than or equal to zero')
  }
  
  N <- -1 # ETC measure using NSRPS is initialized to -1
  L <- length(InputTS)
  
  if (NumBins != 0) {
    SymSeq <- Partition(InputTS, NumBins) # input is a time series of real values
  }  else {
    SymSeq <- InputTS # input is symbolic sequence
  }
  
  # To get rid of zeros in the symbolic sequence, if any
  minY <- min(SymSeq)
  y <- SymSeq - minY
  y <- y + 1
  SymSeq <- y
  
  # The main loop for NSRPS iteration
  N <- 0 # ETC measure 'N'
  Hnew <- ShannonEntropy(SymSeq) # Shannon entropy of the symbolic sequence 
  Len <- length(SymSeq)
  EntLen_array <- Hnew*Len
  
  while ((Hnew > 1e-6) && (Len > 1)) {
    Pair <- FindPair(SymSeq) # find the pair of symbols with maximum frequency
    SymSeqNew <- Substitute(SymSeq, Pair)  # substitute the pair with a new symbol
    Hnew <- ShannonEntropy(SymSeqNew) # Shannon entropy of the new sequence 
    Len <- length(SymSeqNew) 
    EntLen_array <- c(EntLen_array, Hnew*Len)
    N <- N + 1 # ETC measure incremented by 1
    SymSeq <- SymSeqNew
    rm(SymSeqNew)
  }
  
  return(list(N, EntLen_array))
}



# ------------------------------------------------------------------------------
# Subroutine functions
# ------------------------------------------------------------------------------
Partition <- function (InputTS, NumBins) {
  
  # Coverts an input time series into a symbolic sequence
  
  x <- InputTS
  SymSeq <- x
  L <- length(InputTS)
  Range <- max(x) + 1e-6-min(x)
  Delta <- Range / NumBins
  x <- x - min(x)
  
  for (i in 1:L) {
    SymSeq[i] <- floor(x[i] / Delta)
  }
  
  SymSeq <- SymSeq + 1
  
  return(SymSeq)
}


# ------------------------------------------------------------------------------
FindPair <- function (SymSeq) {
  
  # Computes the pair with maximum frequency in the symbolic seq.
  
  Pair <- c()
  Alphabet <- unique(SymSeq)
  M <- max(Alphabet)
  Count_Array <- matrix(0, M, M)
  L <- length(SymSeq)
  indx <- 1
  
  while (indx < L) {
    a <- SymSeq[indx]
    b <- SymSeq[indx+1]
    Count_Array[a, b] <- Count_Array[a, b] + 1
    
    if (a == b) {
      if (indx < (L-1)) {
        if (SymSeq[indx+2] == a) {
          indx <- indx + 1
        }
      }
    }
    indx <- indx + 1
  }
  
  m <- max(Count_Array)
  indx <- which.max(Count_Array)
  Pair1 <- ((indx-1) %% M) + 1
  Pair2 <- floor((indx-1) / M) + 1
  Pair <- c(Pair1, Pair2)
  
  return(Pair)
}



# ------------------------------------------------------------------------------
Substitute <- function (SymSeq, Pair) {
  
  # Pair substitution step of NSRPS
  
  SymSeqNew <- c()
  L <- length(SymSeq)
  Alphabet <- unique(SymSeq)
  M <- length(Alphabet)
  I <- max(Alphabet)
  RepSym <- I + 1 # New Symbol
  indx <- 1
  
  while (indx < L) {
    a <- SymSeq[indx]
    b <- SymSeq[indx + 1]
    
    if (a==Pair[1] && b==Pair[2]) {
      SymSeqNew <- c(SymSeqNew, RepSym)
      indx <- indx + 1
    } else {
      SymSeqNew <- c(SymSeqNew, a)
    }
      
    indx <- indx + 1
    if (indx == L) {
      SymSeqNew <- c(SymSeqNew, SymSeq[length(SymSeq)])
    }
  }
  
  return(SymSeqNew)
}



# ------------------------------------------------------------------------------
ShannonEntropy <- function (SymSeq) {
  
  # Computes the Shannon Entropy of a symbolic sequence

  # To make enteries in the symbolic sequence stictly >0
  minY <- min(SymSeq)
  y <- SymSeq - minY
  y <- y + 1
  
  y <- SymSeq
  H <- 0
  L <- length(y)
  Num <- max(y)
  prob <- matrix(0, 1, Num)
  
  for (i in 1:L) {
    prob[y[i]] <- prob[y[i]] + 1
  }
    prob <- prob / L;
    
    for (i in 1:Num) {
      if (prob[i] != 0) {
        H <- H - prob[i]*log2(prob[i])
      }
    }
    
    return(H)
}

