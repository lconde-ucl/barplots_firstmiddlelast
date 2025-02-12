library(dplyr)
library(ggplot2)
library(stringr)



for (dir in c("human/bigwigs","human/bigwigs_subtract","human/bigwigs_nopseudo", "mouse/bigwigs","mouse/bigwigs_subtract", "mouse/bigwigs_nopseudo")){

	means <- data.frame()

	files<-list.files(dir,pattern=".*input.*txt")   
	for (file in files){
	
		ids = str_split(file, "\\.", simplify=T)
		names<-ids[1]
		name2 = str_split_fixed(names, "_", n=3)
		name<-name2[1]
		tag<-name2[2]
		method<-name2[3]
	
		genome<-ids[4]
		ids2 = str_split_fixed(ids[5], "_", n=2)
		nomalisation<-ids2[1]
		region<-ids2[2]
	
		# mean0 - average over bases with non-covered bases counting as zeroes
   		# mean - average over just covered bases
   		data <- read.table(paste0(dir, "/",file), header = FALSE)
		colnames(data) <- c("name","size","covered","sum","mean0","mean")
		
		# Remove outliers
		# NOTE THIS DOES NOT REMOVE ANYTHINK ON THE LOWER END, AS MOST OF THE TIME THE SIGNAL IS 0
		lower <- quantile(data$mean, 0.01, na.rm=T)
		upper <- quantile(data$mean, 0.99, na.rm=T)	
		
		print(file)
		print("data: ")
		print(dim(data))

		filtered_data <- filter(data, mean >= lower)
		print("data after bottom filtering: ")
		print(dim(filtered_data))

		filtered_data <- filter(filtered_data, mean  <= upper)
		print("data after top filtering: ")
		print(dim(filtered_data))

		# Calculate average signal
		average_signal <- mean(filtered_data$mean)
		
		means <- rbind(means, cbind(name, tag, method, genome, nomalisation, region, average_signal))
	}
	
	means$average_signal<-as.numeric(means$average_signal)
	
	colnames(means)<-c("Sample","Tag", "Method", "Genome","Normalisation","Region","Average_signal")
	
	dir.create(paste0(dir, "_plots"), showWarnings = FALSE)
	write.table(means, paste0(dir, "_plots/average_signal.txt"), sep="\t", row.names=F, col.names=T, quote=F)
	
}
