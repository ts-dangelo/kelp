args <- commandArgs(trailingOnly = TRUE)

indir = args[1]
mdat = args[2]
out = args[3]


if (length(args) <= 1) {
  stop("Input arguments in this order: 1 - directory containing count tables, 2 - metadata file, 3 - absolute path of directory for output", call.=FALSE)
}


library("phyloseq")
library("ANCOMBC")


tables <- list.files(indir, pattern="*.txt", full.names=TRUE)
dir.create(out)
meta <- read.csv(mdat, row.names = 1)
pm <- sample_data(meta)
outliers <- c("Late_Ram_Debbie_S67", "Late_Ram_Otis_S66", "Late_Little_Drisko_BN_S96") # two outlier samples and the blank


for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(benthic)
  	tse$Phase_State = factor(tse$Season, levels = c("May", "August"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Season", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "unpruned_benthic_season_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	benthic_filt = filter_taxa(benthic, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(benthic_filt)
  	tse$Phase_State = factor(tse$Season, levels = c("May", "August"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Season", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "_pruned10_benthic_season_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	may = subset_samples(ps_pr, Season=="May")
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(may)
  	tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "unpruned_benthic_may_phase_state_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	may = subset_samples(ps_pr, Season=="May")
	may_filt = filter_taxa(may, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(may_filt)
  	tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "_pruned10_benthic_may_phase_state_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	August = subset_samples(ps_pr, Season=="August")
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(August)
  	tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "unpruned_benthic_August_phase_state_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	August = subset_samples(ps_pr, Season=="August")
	August_filt = filter_taxa(August, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(August_filt)
  	tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "_pruned10_benthic_August_phase_state_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  	data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  	print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	Midwater = subset_samples(ps_pr, Position=="Midwater")
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(Midwater)
  	tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "_unpruned_Midwater_phase_state_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  	data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  	print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	Midwater = subset_samples(ps_pr, Position=="Midwater")
	Midwater_filt = filter_taxa(Midwater, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(Midwater_filt)
  	tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "_pruned10_Midwater_phase_state_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  	data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  	print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	Midwater = subset_samples(ps_pr, Position=="Midwater")
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(Midwater)
  	tse$Season = factor(tse$Season, levels = c("May", "August"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Season", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "_unpruned_Midwater_Season_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  	data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  	print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	Midwater = subset_samples(ps_pr, Position=="Midwater")
	Midwater_filt = filter_taxa(Midwater, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(Midwater_filt)
  	tse$Season = factor(tse$Season, levels = c("May", "August"))
  
  	set.seed(123)
  	output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Season", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  	res_out = paste0(args[3], "/", data_type, "_pruned10_Midwater_Season_ancom.txt")
  	print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(benthic)
  tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  set.seed(123)
  output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  res_out = paste0(args[3], "/", data_type, "_unpruned_benthic_phase_state_ancom.txt")
  print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}

for (file in tables){
  
  data_type = sapply(strsplit(basename(file),split = "_"), `[`, 1)
  print(data_type)
	
	table <- read.table(file)
	phy_tab <- otu_table(table, taxa_are_rows = FALSE)
	ps_ob <- merge_phyloseq(phy_tab, pm)
	ps_pr <- prune_samples(!(sample_names(ps_ob ) %in% outliers), ps_ob)
	benthic = subset_samples(ps_pr, Position=="Benthic")
	benthic_filt = filter_taxa(benthic, function(x) sum(x > 0) > (0.1*length(x)), TRUE)
	
	tse = mia::makeTreeSummarizedExperimentFromPhyloseq(benthic_filt)
  tse$Phase_State = factor(tse$Phase_State, levels = c("Kelp", "Turf"))
  
  set.seed(123)
  output = ancombc2(data = tse, assay_name = "counts",
                  fix_formula = "Phase_State", rand_formula = NULL,
                  p_adj_method = "holm", pseudo_sens = TRUE,
                  prv_cut = 0.10, lib_cut = 1000, s0_perc = 0.05, n_cl = 8)
  
  res_out = paste0(args[3], "/", data_type, "_pruned10_benthic_phase_state_ancom.txt")
  print(res_out)
	write.table(output$res, res_out, sep = "\t", col.names=NA)
}
