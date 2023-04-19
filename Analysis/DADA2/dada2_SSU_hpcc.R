library(dada2)

path <- "/mnt/scratch/f0008425/SSU"

set.seed(1390)
CORES = parallel::detectCores()

# Forward and reverse fastq filenames have the format:
fnFs <- sort(list.files(path, pattern = "_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format:
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(fnFs, get.sample.name))
head(sample.names)

# Filter and trim
filtFs <- file.path(path, "filtered", basename(fnFs))
filtRs <- file.path(path, "filtered", basename(fnRs))
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240), maxN = 0, maxEE = c(2, 2), truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = CORES)

# If the number of sequence files were different after filtering
path_edited <- "/mnt/scratch/f0008425/SSU/filtered"

fnFs_edited <- sort(list.files(path_edited, pattern = "_R1.fastq", full.names = TRUE))
fnRs_edited <- sort(list.files(path_edited, pattern = "_R2.fastq", full.names = TRUE))
errF <- learnErrors(fnFs_edited, multithread = CORES)
errR <- learnErrors(fnRs_edited, multithread = CORES)

# Dereplicate identical reads
derepFs <- derepFastq(fnFs_edited, verbose = TRUE)
derepRs <- derepFastq(fnRs_edited, verbose = TRUE)

# Name the derep-class objects by the sample names
sample.names_edited <- unname(sapply(fnFs_edited, get.sample.name))
names(derepFs) <- sample.names_edited
names(derepRs) <- sample.names_edited

# Sample Inference
dadaFs <- dada(derepFs, err = errF, multithread = CORES)
dadaRs <- dada(derepRs, err = errR, multithread = CORES)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

# Construct Sequence Table
seqtab <- makeSequenceTable(mergers)

# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=CORES, verbose=TRUE)

# Track reads through the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# Assign taxanomy
silva.ref <- "/mnt/scratch/f0008425/SSU/silva_nr99_v138.1_train_set.fa.gz"
silvaspecies.ref <- "/mnt/scratch/f0008425/SSU/silva_species_assignment_v138.1.fa.gz"
taxa <- assignTaxonomy(seqtab.nochim, silva.ref, multithread = CORES, tryRC = TRUE)
taxa <- addSpecies(taxa, silvaspecies.ref)

saveRDS(t(seqtab.nochim), paste0(path_edited, "/seq_table.RDS"))
saveRDS(taxa, paste0(path_edited, "/tax_table.RDS"))
