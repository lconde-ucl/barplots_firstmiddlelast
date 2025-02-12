library(dplyr)
library(ggplot2)
library(stringr)



for (dir in c("human/bigwigs", "human/bigwigs_subtract", "human/bigwigs_nopseudo", "mouse/bigwigs", "mouse/bigwigs_subtract", "mouse/bigwigs_nopseudo")){

	means<-read.table(paste0(dir, "_plots/average_signal.txt"), sep="\t", header=T)

	# Reorder the Region
	means$Region <- factor(means$Region, levels = c("upTSS", "first_exon", "first_intron", "middle_exon", "middle_intron", "last_intron", "last_exon","dnTES"))
	
	png(paste0(dir, "_plots/average_signal.png"), width = 600, height = 450)
	p<-ggplot(means, aes(x = Region, y = Average_signal, color = Sample, group = Sample)) + 
	    geom_line() +  # Connect points with lines
            geom_point() +  # Optionally add points to the plot
	    theme_minimal() +
	    labs(title = "Average Profiles by Region (Faceted by Method and Tag)",
	       x = "Gene Regions",
	       y = "Average Signal (CPM)") +
	  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	  facet_wrap(Method ~ Tag)  # Create panels for each combination of Method and Tag
	print(p)
	dev.off()

	for (sample in unique(means$Sample)){
		for (method in unique(means$Method)){
			means_sample<-filter(means, Sample==sample & Method==method)
			png(paste0(dir, "_plots/average_signal_", sample, "_", method, ".png"), width = 600, height = 450)
			p <- ggplot(means_sample, aes(x = Region, y = Average_signal, color = Tag, fill = Tag)) + 
  				geom_bar(stat = "identity", position = "dodge") +  # Use bars instead of lines
  				theme_minimal() +
  				labs(title = paste("Average Profiles by Region -", sample, method),
       				x = "Gene Regions",
       				y = "Average Signal (CPM)") +
  				theme(axis.text.x = element_text(angle = 45, hjust = 1))
			print(p)
			dev.off()
		}
	}

}

