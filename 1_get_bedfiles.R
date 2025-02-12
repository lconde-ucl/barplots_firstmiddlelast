library(dplyr)
library(stringr)


#- Human genome
#canonical_transcript="knownCanonical.txt"
#exons_file="wgEncodeGencodeBasicV40.tsv"
##exons_file="exons.txt" #- for testing
#genome="human/hs"

#- Mouse genome
canonical_transcript="knownCanonical.Mm.txt"
exons_file="wgEncodeGencodeBasicVM23.tsv"
genome="mouse/mm"


dir.create(paste0(genome, "_bedfiles"), recursive=T)

#- get 1 transcript per gene only (remove non canonical chromsomes)
canonical<-read.table(canonical_transcript, sep="\t", header=F)
colnames(canonical)<-c("chr","start","end","id","transcript_id","gene_id")
dim(canonical)
canonical<-canonical[!str_detect(canonical$chr, "_"),]
canonical <- filter(canonical, chr != "chrM")
dim(canonical)


#- get exons in anonical transcripts only
#- some transcripts appear in both chrX and chrY in the exons file (eg  ENST00000627721.3)
#- so search for transcript-chr pairs
transcripts<-paste0(canonical$transcript_id, "-", canonical$chr)
exons<-read.table(exons_file, header=T, sep="\t")
dim(exons)
exons<-filter(exons, paste0(name, "-",chrom)  %in% transcripts)
dim(exons)
dim(unique(exons))


first_exon <- data.frame()
last_exon <- data.frame()
first_intron <- data.frame()
last_intron <- data.frame()
upTSS <- data.frame()
dnTES <- data.frame()
middle_exon <- data.frame()
middle_intron <- data.frame()

for (i in c(1:length(exons$name))){

  d<-exons[i,]

  start_positions <- as.numeric(unlist(strsplit(d$exonStarts, ",")))
  end_positions <- as.numeric(unlist(strsplit(d$exonEnds, ",")))

  if(d$strand == "+"){
    #print(d$name)

    first_exon <- rbind(first_exon, cbind(d$chrom, sprintf("%.0f", start_positions[1]), sprintf("%.0f", end_positions[1]), d$name))
    last_exon <- rbind(last_exon, cbind(d$chrom, sprintf("%.0f", start_positions[d$exonCount]), sprintf("%.0f", end_positions[d$exonCount]), d$name))
 
    upTSS <- rbind(upTSS, cbind(d$chrom, sprintf("%.0f", pmax(d$txStart-1000, 0)), sprintf("%.0f", d$txStart-1), d$name))
    dnTES <- rbind(dnTES, cbind(d$chrom, sprintf("%.0f", d$txEnd+1), sprintf("%.0f", d$txEnd+1000), d$name))
    
    #- if there is 1 exon, this will appear in first_exon and last_exon, but there will be no middle_exon, nor first intron
    #- if there are 2 exons, there will be first_intron (which will be the same as last intron), but still no middle_exon
    #- if there are 3 exons, there will be first and last intron, as well as middle_exons, but still no middle intron
    #- if there are more than 3 exons, there will be middle_exons and middle intron
    if(d$exonCount > 1) {
      first_intron <- rbind(first_intron, cbind(d$chrom, sprintf("%.0f", end_positions[1]+1), sprintf("%.0f", start_positions[2]-1), d$name))
      last_intron <- rbind(last_intron, cbind(d$chrom, sprintf("%.0f", end_positions[d$exonCount-1]+1), sprintf("%.0f", start_positions[d$exonCount]-1), d$name))
    }

  }else{

    #- If negative strand, switch first with last exon
    first_exon <- rbind(first_exon, cbind(d$chrom, sprintf("%.0f", start_positions[d$exonCount]), sprintf("%.0f", end_positions[d$exonCount]), d$name))
    last_exon <- rbind(last_exon, cbind(d$chrom, sprintf("%.0f", start_positions[1]), sprintf("%.0f", end_positions[1]), d$name))

    upTSS <- rbind(upTSS, cbind(d$chrom, sprintf("%.0f", d$txEnd+1), sprintf("%.0f", d$txEnd+1000), d$name))
    dnTES <- rbind(dnTES, cbind(d$chrom, sprintf("%.0f", pmax(d$txStart-1000, 0)), sprintf("%.0f", d$txStart-1), d$name))

    if(d$exonCount > 1) {
      #- If negative strand, switch first with last intron
      first_intron <- rbind(first_intron, cbind(d$chrom, sprintf("%.0f", end_positions[d$exonCount-1]+1), sprintf("%.0f", start_positions[d$exonCount]-1), d$name))
      last_intron <- rbind(last_intron, cbind(d$chrom, sprintf("%.0f", end_positions[1]+1), sprintf("%.0f", start_positions[2]-1), d$name))
    }
  }
  if(d$exonCount > 2) {
    for (j in c(2:(length(start_positions)-1))){
      middle_exon <- rbind(middle_exon, cbind(d$chrom, sprintf("%.0f", start_positions[j]), sprintf("%.0f", end_positions[j]), paste0(d$name, "-", j)))
    }
  }
  if(d$exonCount > 3) {
    for (k in c(2:(length(start_positions)-2))){
      middle_intron <- rbind(middle_intron, cbind(d$chrom, sprintf("%.0f", end_positions[k]+1), sprintf("%.0f", start_positions[k+1]-1), paste0(d$name, "-", k)))
    }
  }
}


write.table(first_exon, paste0(genome, "_bedfiles/first_exon.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(last_exon, paste0(genome, "_bedfiles/last_exon.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(first_intron, paste0(genome, "_bedfiles/first_intron.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(last_intron, paste0(genome, "_bedfiles/last_intron.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(upTSS, paste0(genome, "_bedfiles/upTSS.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(dnTES, paste0(genome, "_bedfiles/dnTES.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(middle_exon, paste0(genome, "_bedfiles/middle_exon.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(middle_intron, paste0(genome, "_bedfiles/middle_intron.bed"), sep="\t", row.names=F, col.names=F, quote=F)


