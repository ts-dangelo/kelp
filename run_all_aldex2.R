args <- commandArgs(trailingOnly = TRUE)

indir = args[1]
mdat = args[2]
out = args[3]


if (length(args) <= 1) {
  stop("Input arguments in this order: 1 - directory containing count tables, 2 - metadata file, 3 - absolute path of directory for output", call.=FALSE)
}


library("phyloseq")
library("compositions")
library("ALDEx2")
#remove.packages("tidyverse");install.packages("tidyverse")


tables <- list.files(indir, pattern="*.txt", full.names=TRUE)
dir.create(out)
meta <- read.csv(mdat, row.names = 1)
pm <- sample_data(meta)
outliers <- c("Late_Ram_Debbie_S67", "Late_Ram_Otis_S66", "Late_Little_Drisko_BN_S96") # two outlier samples and the blank

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	
	
	ps_T <- t(otu_table(ps_pr))
	pr_mdat <- sample_data(ps_pr)
	conds <- as.character(pr_mdat$Position)
	
	res_out = paste0(args[3], "/", data_type, "_unpruned_position_aldex2.txt")
	print(res_out)
	
	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")
    
	write.table(results, res_out, sep = "\t", col.names=NA)

}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	ps_prev = filter_taxa(ps_pr, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	ps_T <- t(otu_table(ps_prev))
	pr_mdat <- sample_data(ps_prev)
	conds <- as.character(pr_mdat$Position)
	
	res_out = paste0(args[3], "/", data_type, "_pruned10_position_aldex2.txt")
	print(res_out)
	
	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	
	
	ps_T <- t(otu_table(benthic))
	pr_mdat <- sample_data(benthic)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type, "_unpruned_benthic_phase_state_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	benthic_prev = filter_taxa(benthic, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	ps_T <- t(otu_table(benthic_prev))
	pr_mdat <- sample_data(benthic_prev)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type, "_pruned10_benthic_phase_state_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	
	
	ps_T <- t(otu_table(benthic))
	pr_mdat <- sample_data(benthic)
	conds <- as.character(pr_mdat$Season)
	
	res_out = paste0(args[3], "/", data_type, "_unpruned_benthic_season_aldex2.txt")
	print(res_out)
	
	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	benthic_prev = filter_taxa(benthic, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	ps_T <- t(otu_table(benthic_prev))
	pr_mdat <- sample_data(benthic_prev)
	conds <- as.character(pr_mdat$Season)
	
	res_out = paste0(args[3], "/", data_type, "_pruned10_benthic_season_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	august = subset_samples(ps_pr, Season=="August")
	
	
	ps_T <- t(otu_table(august))
	pr_mdat <- sample_data(august)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type, "_unpruned_august_benthic_phase_state_aldex2.txt")
	print(res_out)
	
	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	august = subset_samples(ps_pr, Season=="August")
	august_prev = filter_taxa(august, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	ps_T <- t(otu_table(august_prev))
	pr_mdat <- sample_data(august_prev)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type, "_pruned10_august_benthic_phase_state_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	may = subset_samples(ps_pr, Season=="May")
	
	
	ps_T <- t(otu_table(may))
	pr_mdat <- sample_data(may)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type,"_unpruned_may_benthic_phase_state_aldex2.txt")
	print(res_out)
	
	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	may = subset_samples(ps_pr, Season=="May")
	may_prev = filter_taxa(may, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	ps_T <- t(otu_table(may_prev))
	pr_mdat <- sample_data(may_prev)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type,"_pruned10_may_benthic_phase_state_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	midwater = subset_samples(ps_pr, Position=="Midwater")
	
	ps_T <- t(otu_table(midwater))
	pr_mdat <- sample_data(midwater)
	conds <- as.character(pr_mdat$Season)
	
	res_out = paste0(args[3], "/", data_type,"_unpruned_midwater_season_aldex2.txt")
	print(res_out)


	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	midwater = subset_samples(ps_pr, Position=="Midwater")
	midwater_prev = filter_taxa(midwater, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	ps_T <- t(otu_table(midwater_prev))
	pr_mdat <- sample_data(midwater_prev)
	conds <- as.character(pr_mdat$Season)
	
	res_out = paste0(args[3], "/", data_type,"_pruned10_midwater_season_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	midwater = subset_samples(ps_pr, Position=="Midwater")
	
	ps_T <- t(otu_table(midwater))
	pr_mdat <- sample_data(midwater)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type,"_unpruned_midwater_phase_state_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	midwater = subset_samples(ps_pr, Position=="Midwater")
	midwater_prev = filter_taxa(midwater, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	ps_T <- t(otu_table(midwater_prev))
	pr_mdat <- sample_data(midwater_prev)
	conds <- as.character(pr_mdat$Phase_State)
	
	res_out = paste0(args[3], "/", data_type,"_pruned10_midwater_phase_state_aldex2.txt")
	print(res_out)

	results <- aldex(reads=ps_T, conditions = conds, mc.samples = 128, test="t", effect=TRUE,
                 include.sample.summary = FALSE, verbose=T, denom="all")

	write.table(results, res_out, sep = "\t", col.names=NA)
}

