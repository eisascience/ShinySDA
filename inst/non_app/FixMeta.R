
base.path = "/Volumes/Maggie/Work/OHSU/Tejpal/TejP_hu_Autoimm/data/sda/primeseq/"
list.files(base.path, pattern = ".rds")

SDA_MetaDF <- readRDS(paste0(base.path, "sda.566999/566999_MetaOldDF.rds"))
rownames(SDA_MetaDF) = SDA_MetaDF$cellbarcode

MetaDF_full <- readRDS("/Volumes/Maggie/Work/OHSU/Tejpal/TejP_hu_Autoimm/data/FullCombo_MetaDF.rds")





SDA_MetaDF = SDA_MetaDF[,!grepl("celltypist", colnames(SDA_MetaDF))]
SDA_MetaDF = SDA_MetaDF[,!grepl("UCell", colnames(SDA_MetaDF))]
SDA_MetaDF = SDA_MetaDF[,!grepl("is.pure", colnames(SDA_MetaDF))]


MetaOfInterest = c("scGateConsensus", "BarcodePrefix", "immgen.label", "monaco.label",
                   "SingleRConsensus", "cDNA_ID", "SampleDate", "SubjectId", "barcode",
                   "Tissue", "Population", "BatchId", "Age", "Sex", "Race", "Final.Dx", 
                   "Dx", "Addl.Dx", "B27pos", "Systemic.Activity", "Eye.Activity", "CellViability", "Pheno1")



commonFeats <- intersect(colnames(SDA_MetaDF), colnames(MetaDF_full) )


SDA_MetaDF_new <- cbind(SDA_MetaDF[,!(colnames(SDA_MetaDF) %in% commonFeats)], MetaDF_full[rownames(SDA_MetaDF), ])

SDA_MetaDF_new$SubjectId

saveRDS(SDA_MetaDF_new, paste0(base.path, "sda.566999/566999_MetaDF.rds"))
