# File initially created on Monday, September 21, 2020

# This file contains a collection of functions used to transform a Sniffles VCF file
# into a format that will be usable by bcftools norm in order to then merge different
# vcf files for genotyping with programs such as BayesTyper/vg/Paragraph

# As of now, there seems to be a problem with Sniffles as the reference alleles reported
# do not match the coordinates reported for the variant. Here, I will consider that the
# coordinates that are reported are right and will update the reference and alternate alleles
# accordingly.

# The main function (fix_sniffles) takes care of input/output and extracts relevant
# information from the vcf file that can be passed to other functions that compute the
# information that will be used to update the fields of the VCF. There will be one
# such function for each type of SV considered (DEL, INS, INV, DUP).

# WARNING: The two following lines are no longer true in the current version of the function
# FORMAT and genotype fields, as well as unnecessary INFO fields will be stripped from
# the vcf such that only information that is necessary for downstream genotyping is used.

# input_vcf: a character string with the name of the input vcf file
# output_vcf: a character string with the name of the output vcf file (WARNING: will overwrite file)
# refgenome: a character string with the location of the file corresponding to the reference genome
#            in FASTA format. A .fai index matching the file must also exist.
fix_sniffles <- function(input_vcf, output_vcf, refgenome = NULL) {

	# Checking if the reference genome has been supplied
	if(is.null(refgenome) || !file.exists(refgenome)) {
		stop("Reference genome file is not supplied or does not exist")
	}

	# Checking if the reference genome is indexed
	if(!file.exists(paste0(refgenome, ".fai"))) {
		stop("Reference genome .fai index does not exist")
	}

	# Opening a connection to read the vcf file
	input_con <- file(input_vcf, open = "rt")
	on.exit(close(input_con), add = TRUE)

	# Opening another connection to the output file
	output_con <- file(output_vcf, open = "wt")
	on.exit(close(output_con), add = TRUE)

	# Lines are then used one by one to output the header
	while(grepl("^#", cur_line <- scan(input_con, what = character(), sep = "\n", n = 1, quiet = TRUE))) {
		cat(paste0(cur_line, "\n"), file = output_con)
	}

	# Reading the vcf from file wnad ignoring header lines
	vcf <- read.table(input_vcf, comment.char = "#", stringsAsFactors = FALSE)

	# Creating indices vectors for each SV type
	dels <- which(grepl("SVTYPE=DEL", vcf[[8]]))
	ins  <- which(grepl("SVTYPE=INS", vcf[[8]]))
	dups <- which(grepl("SVTYPE=DUP", vcf[[8]]))
	invs <- which(grepl("SVTYPE=INV", vcf[[8]]))
	
	# Extracting some useful information for each variant
	chrs   <- vcf[[1]]
	starts <- vcf[[2]]
	altseq <- vcf[[5]]
	svlen  <- as.numeric(sub(".*SVLEN=(-?[0-9]+);.*", "\\1", vcf[[8]]))
	#svtype <- sub(".*(SVTYPE=[A-Z]+);.*", "\\1", vcf[[8]])
	#imprecise <- ifelse(grepl("IMPRECISE", vcf[[8]]), "IMPRECISE;", "")
	#realigned <- ifelse(grepl("REALIGNED", vcf[[8]]), "REALIGNED;", "")

	# Computing the replacement information for each variant
	del_info <- del_process(chrs[dels], starts[dels], abs(svlen[dels]), refgenome = refgenome)
	ins_info <- ins_process(chrs[ins],  starts[ins],  altseq[ins],      refgenome = refgenome)
	dup_info <- dup_process(chrs[dups], starts[dups], abs(svlen[dups]), refgenome = refgenome)
	#inv_info <- inv_process(chrs[invs], starts[invs], abs(svlen[invs]), refgenome = refgenome)

	# Assigning the results to the right columns of the vcf file
	vcf[[2]][dels] <- del_info$pos
	vcf[[4]][dels] <- del_info$ref
	vcf[[5]][dels] <- del_info$alt

	vcf[[2]][ins]  <- ins_info$pos
	vcf[[4]][ins]  <- ins_info$ref
	vcf[[5]][ins]  <- ins_info$alt

	vcf[[2]][dups] <- dup_info$pos
	vcf[[4]][dups] <- dup_info$ref
	vcf[[5]][dups] <- dup_info$alt

	#vcf[[4]][invs] <- inv_info$ref
	#vcf[[5]][invs] <- inv_info$alt

	# Reformatting the INFO column with only the information we need
	# vcf[[8]] <- paste0(imprecise, realigned, svtype, ";SVLEN=", as.character(svlen))

	# Writing the data.frame to the output file
	write.table(vcf, file = output_con, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

	# Writing
	return(invisible(NULL))
}

del_process <- function(chr, start, width, refgenome) {

	# The deletion start needs to be before the actual deletion
	pos <- ifelse(start == 1, 1, start - 1)

	# The last deleted nucleotide is width units further
	end <- pos + width

	# Creating a GRanges object for querying the ref sequence
	ref_range <- GenomicRanges::GRanges(seqnames = chr, 
					    ranges = IRanges::IRanges(start = pos, end = end))

	# Querying the reference sequence
	ref <- Rsamtools::scanFa(refgenome, ref_range)

	# Creating a GRanges object for the alt sequence
	alt_range <- GenomicRanges::GRanges(seqnames = chr,
					    ranges = IRanges::IRanges(start = pos, end = pos))

	# Querying the alt sequence
	alt <- Rsamtools::scanFa(refgenome, alt_range)

	# Returning the formatted information
	list(pos = pos, ref = unname(as.character(ref)), alt = unname(as.character(alt)))
}

ins_process <- function(chr, start, alt_seq, refgenome) {

	# The start of the insertion must be offset by one
	pos <- ifelse(start == 1, 1, start - 1)

	# Creating a GRanges object for querying the reference sequence
	ref_range <- GenomicRanges::GRanges(seqnames = chr,
					    ranges = IRanges::IRanges(start = pos, end = pos))

	# Querying the reference nucleotide
	ref <- Rsamtools::scanFa(refgenome, ref_range)
	ref <- unname(as.character(ref))

	# The alternate sequence is the reference sequence to which the alt_seq from Sniffles is added
	alt <- paste0(ref, alt_seq)

	# Returning the formatted information
	list(pos = pos, ref = ref, alt = alt)
}

dup_process <- function(chr, start, width, refgenome) {

	# The start of the duplication must be offset by one
	pos <- ifelse(start == 1, 1, start - 1)

	# The last duplicated nucleotide is width units further
	end <- pos + width

	# Creating a GRanges object for querying the alt sequence
	alt_range <- GenomicRanges::GRanges(seqnames = chr,
					    ranges = IRanges::IRanges(start = pos, end = end))

	# Querying the alt sequence
	alt <- Rsamtools::scanFa(refgenome, alt_range)

	# Creating a GRanges object for querying the reference sequence
	ref_range <- GenomicRanges::GRanges(seqnames = chr,
					    ranges = IRanges::IRanges(start = pos, end = pos))

	# Querying the reference sequence
	ref <- Rsamtools::scanFa(refgenome, ref_range)

	# Returning the formatted information
	list(pos = pos, ref = as.character(unname(ref)), alt = as.character(unname(alt)))
}

inv_process <- function(chr, start, width, refgenome) {

	# In the inversion case, the variation truly starts at "start" so we need not update the position
	# However, the end position will be start + width -1
	end <- start + width - 1

	# Creating a GRanges object to query the reference sequence
	ref_range <- GenomicRanges::GRanges(seqnames = chr,
					    ranges = IRanges::IRanges(start = start, end = end))

	# Querying the reference sequence
	ref <- Rsamtools::scanFa(refgenome, ref_range)

	# The alt sequence is simply the reverse complement
	alt <- sapply(unname(as.character(ref)), revcomp, USE.NAMES = FALSE)

	list(pos = start, ref = as.character(unname(ref)), alt = alt)
}

revcomp <- function(sequence) {

	# Checking that only one sequence is provided
	stopifnot(length(sequence) == 1)

	# The lookup table that will be used for replacement
	rep_table <- c("A" = "T",
		       "T" = "A",
		       "G" = "C",
		       "C" = "G",
		       "N" = "N")

	# Splitting the sequence into its constituent nucleotides
	sequence <- strsplit(sequence, "")[[1]]

	# Replacing the nucleotides by their complement
	sequence <- rep_table[sequence]

	# Returning the inverted sequence
	paste0(rev(sequence), collapse = "")
}

