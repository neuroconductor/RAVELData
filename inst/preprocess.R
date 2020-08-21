library(kirby21.scan.1)
library(pbapply)
library(ANTsR)
library(extrantsr)
library(fslr)
library(RAVEL)
library(EveTemplate)

# Getting EveTemplate data first:
template   <- antsImageRead(EveTemplate::getEvePath(),3) # Eve template:
brain_mask <- antsImageRead(EveTemplate::getEvePath("Brain_Mask"), 3)

# Getting KIRBY21 data:
kirbyDir  <- system.file(package = "kirby21.scan.1")
files <- list.files(kirbyDir,
	recursive=TRUE,
	pattern="MPRAGE.nii.gz",
	full.names=TRUE
)[1:4]

# Reading the files:
scans <- pblapply(files, function(x) antsImageRead(x, 3))
# N4 Correction:
scans_n4 <- pblapply(scans, n4BiasFieldCorrection)
# Non-linear registration:
scans_n4_reg <- pblapply(scans_n4, function(x){
	warp <- antsRegistration(fixed = template, moving = x, typeofTransform = "SyN")
	temp <- antsImageClone(warp$warpedmovout)
	temp
})
# Brain extraction:
scans_n4_reg_brain <- pblapply(scans_n4_reg, function(x){
	x[brain_mask==0] <- 0
	x
})
# Segmentation with FSL FAST:
scans_n4_reg_brain <- pblapply(scans_n4_reg_brain, ants2oro)
segs  <- pblapply(scans_n4_reg_brain, function(x){
	fast(file=x, outfile="try", opts= "-t 1 -n 3 --nobias", retimg=TRUE, verbose=FALSE)	
})
system("rm try*")
segs_csf <- pblapply(segs, function(x){
	x[x!=1] <- 0
	x
})
# Writing out the processed images:
pblapply(1:4, function(j){
	prefix <- "extdata/scan"
	antsImageWrite(scans[[j]], paste0(prefix, j, ".nii.gz"))
	antsImageWrite(scans_n4_reg_brain[[j]], paste0(prefix, j, "_processed.nii.gz"))
	writeNIfTI(segs[[j]], paste0(prefix, j, "_seg"))
	writeNIfTI(segs_csf[[j]], paste0(prefix, j, "_csf_mask"))
})

#Let's create the mask intersect:
masks <- list.files("extdata", full.names=TRUE, pattern="*mask*.nii*")
mask  <- maskIntersect(masks, output.file="extdata/csf_mask_intersection.nii.gz")


# Alternative segmentation with atropos from ANTs:
#scans_n4_reg_brain <- pblapply(scans_n4_reg_brain, oro2ants)
#segs  <- pblapply(scans_n4_reg_brain, otropos)
#segs_csf <- pblapply(segs, function(x){
#  x[x!=1] <- 0
#  x
#})
#pblapply(1:4, function(j){writeNIfTI(segs[[j]], paste0("extdata/scan",j,"_seg"))})
#pblapply(1:4, function(j){writeNIfTI(segs_csf[[j]], paste0("extdata/scan",j,"_csf_mask"))})


#pdf("try.pdf", width=5, height=5)
#multi_overlay(scans_n4_reg_brain)
#dev.off()
# multi_overlay(scans_n4_reg_brain, 
#   par.opts=list(mfrow=c(2,2),bg = "black",oma = c(0, 0, 0, 0), mar = rep(0, 4)))
# dev.off()

# files <- list.files("../data/", pattern="processed.nii.gz", full.names=TRUE)
# scans_n4_reg_brain <- lapply(files, readNIfTI)


