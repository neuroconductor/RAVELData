library(kirby21.scan.1)
library(pbapply)
library("ANTsR", lib.loc="/dexter/disk2/smart/programs/R_Packages")
library(fslr)




# Written by John Muschelli
ants2oro <- function(img, 
                     reorient = FALSE){
  if ( is.antsImage(img) | is.character(img) ) {
    fname = tempants(img)
    img = readNIfTI(fname, reorient = reorient)
    return(img)
  }
  if ( is.nifti(img) ) {
    return(img)
  }
  stop("img not class nifti or antsImage")
  return(NULL)
}

# Written by John Muschelli
tempants <- function(x, # object of class \code{antsImage}
                     gzipped = TRUE # logical if the tempfile should be gzipped or not
){
  if (inherits(x, "character")) {
    return(x)
  } else {
    if (inherits(x, "antsImage")) {
      ext = ".nii"
      if (gzipped) ext = ".nii.gz"
      tfile = paste0(tempfile(), ext)
      antsImageWrite(x, tfile)
      return(tfile)
    }  else {
      stop("x has unknown class - not char or nifti")
    }
  }
  return(FALSE)
}





kirbyDir  <- system.file(package = "kirby21.scan.1")
files     <- list.files(kirbyDir, recursive=TRUE, pattern="MPRAGE.nii.gz", full.names=TRUE)[1:4]

# Reading the files:
scans <- pblapply(files, function(x) antsImageRead(x, 3))

# N4 Correction:
scans_n4 <- pblapply(scans, n4BiasFieldCorrection)

# Non-linear registration:
template  <- antsImageRead("../data/JHU_MNI_SS_T1.nii.gz",3) # Eve template:
scans_n4_reg <- pblapply(scans_n4, function(x){
	warp <- antsRegistration(fixed = template, moving = x, typeofTransform = "SyN")
	temp <- antsImageClone(warp$warpedmovout)
	temp
})

# Brain extraction:
brain_mask <- antsImageRead("../data/JHU_MNI_SS_T1_Brain_Mask.nii.gz", 3)
scans_n4_reg_brain <- pblapply(scans_n4_reg, function(x){
	x[brain_mask==0] <- 0
	x
})

# Writing out the processed images:
pblapply(1:4, function(j){
	antsImageWrite(scans_n4_reg_brain[[j]], 
		paste0("../data/scan",j,"_processed.nii.gz"))
})


scans_n4_reg_brain <- pblapply(scans_n4_reg_brain, ants2oro)

# Segmentation:
segs  <- pblapply(scans_n4_reg_brain, function(x){
	fast(file=x, outfile="try", opts= "-t 1 -n 3 --nobias", retimg=TRUE, verbose=FALSE)	
})
system("rm try*")

pblapply(1:4, function(j){
	writeNIfTI(segs[[j]], 
		paste0("../data/scan",j,"_seg.nii.gz"))
})

segs_csf <- pblapply(segs, function(x){
	x[x!=1] <- 0
	x
})


pblapply(1:4, function(j){
	writeNIfTI(segs_csf[[j]], 
		paste0("../data/scan",j,"_csf_mask.nii.gz"))
})












