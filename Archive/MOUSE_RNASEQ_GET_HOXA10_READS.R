library(Rsubread)

bams = list.files(pattern=".bam")

res1=featureCounts(bams[1], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
saveRDS(res1, file=paste(bams[1], "counts.rds", sep=""))

res2=featureCounts(bams[2], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
saveRDS(res2, file=paste(bams[2], "counts.rds", sep=""))

res3=featureCounts(bams[3], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
saveRDS(res3, file=paste(bams[3], "counts.rds", sep=""))

res4=featureCounts(bams[4], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
saveRDS(res4, file=paste(bams[4], "counts.rds", sep=""))

res5=featureCounts(bams[5], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
saveRDS(res5, file=paste(bams[5], "counts.rds", sep=""))

res6=featureCounts(bams[6], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
saveRDS(res6, file=paste(bams[6], "counts.rds", sep=""))

#res7=featureCounts(bams[7], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
#saveRDS(res7, file=paste(bams[1], "counts.rds", sep=""))

res8=featureCounts(bams[8], isGTFAnnotationFile = TRUE, annot.ext = "hoxa_region.gtf", GTF.featureType = "transcript", isPairedEnd=TRUE)
saveRDS(res8, file=paste(bams[8], "counts.rds", sep=""))

