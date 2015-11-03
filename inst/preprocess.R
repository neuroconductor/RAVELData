library(kirby21.scan.1)
library(pbapply)
library("ANTsR", lib.loc="/dexter/disk2/smart/programs/R_Packages")

kirbyDir  <- system.file(package = "kirby21.scan.1")
files     <- list.files(kirbyDir, recursive=TRUE, pattern="MPRAGE.nii.gz", full.names=TRUE)[1:4]

# Reading the files:
scans <- pblapply(files, function(x) antsImageRead(x, 3))

# N4 Correction:
scans_n4 <- pblapply(scans, n4BiasFieldCorrection)

# Non-linear registration:
template  <- antsImageRead("../atlases/JHU_MNI_SS_T1.nii.gz",3) # Eve template:
scans_n4_reg <- pblapply(scans_n4, function(x){
	warp <- antsRegistration(fixed = template, moving = x, typeofTransform = "SyN")
	temp <- antsImageClone(warp$warpedmovout)
	temp
})

# Brain extraction:
brain_mask <- antsImageRead("../atlases/JHU_MNI_SS_T1_Brain_Mask.nii.gz", 3)
scans_n4_reg_brain <- pblapply(scans_n4_reg, function(x){
	x[brain_mask==0] <- 0
	x
})

# Writing out the processed images:
pblapply(1:4, function(j){
	antsImageWrite(scans_n4_reg_brain[[j]], 
		paste0("../processed/scan",j,"_processed.nii.gz"))
})









