# Build the GMT file from my curated gene set workbook
#
# The workbook (MRS_Metabolite_GeneSets_v4.xlsx) has one row per gene set:
#   Metabolite | Reason | Gene Set (MSigDB) | Gene Count | Link | Gene 1 ... Gene n
# The metabolite name is only written on the first row of each block.
#
# GSEA needs a .gmt file instead: one line per gene set, tab separated,
#   name <tab> description <tab> gene1 <tab> gene2 ...
# Every gene set name ends in a metabolite suffix (_CHO, _MI, _NAA ...), which is
# how the later scripts know which metabolite a set belongs to.
#
# Run from the repository root:  Rscript scripts/01_build_gene_set_gmt.R

library(readxl)

xlsx_file <- "data/gene_sets/MRS_Metabolite_GeneSets_v4.xlsx"
gmt_file  <- "data/gene_sets/MRS_Metabolite_GeneSets_suffix.gmt"

tab <- as.data.frame(read_excel(xlsx_file, sheet = 1, col_types = "text"))
tab <- tab[!is.na(tab[[3]]), ]            # keep rows that actually hold a gene set
cat("Gene sets in workbook:", nrow(tab), "\n")

set_names <- trimws(tab[[3]])
links     <- ifelse(is.na(tab[[5]]), "na", tab[[5]])

# gene set names have to be unique or GSEA silently drops the repeats
dups <- set_names[duplicated(set_names)]
if (length(dups) > 0) stop("Duplicate gene set names: ", paste(dups, collapse = ", "))

lines <- character(nrow(tab))
n_genes <- integer(nrow(tab))
for (i in seq_len(nrow(tab))) {
  genes <- trimws(unlist(tab[i, 6:ncol(tab)]))
  genes <- unique(genes[!is.na(genes) & genes != ""])
  n_genes[i] <- length(genes)
  lines[i] <- paste(c(set_names[i], links[i], genes), collapse = "\t")
}
writeLines(lines, gmt_file)

# check the Gene Count column still matches what is actually listed
declared <- as.integer(tab[[4]])
if (any(declared != n_genes)) {
  cat("Gene Count column disagrees for:\n")
  print(set_names[declared != n_genes])
}

# Excel likes turning gene symbols into dates (SEPT2 -> 2-Sep), so check for that
all_genes <- unique(unlist(strsplit(sub("^[^\t]*\t[^\t]*\t", "", lines), "\t")))
date_like <- grep("^[0-9]{1,2}-[A-Za-z]{3}$", all_genes, value = TRUE)
if (length(date_like) > 0) cat("Date-converted symbols found:", date_like, "\n")

suffix <- sub(".*_", "", set_names)
cat("\nGene sets per metabolite suffix:\n")
print(table(suffix))
cat(sprintf("\nWrote %d gene sets (%d unique genes, set sizes %d-%d) to %s\n",
            length(lines), length(all_genes), min(n_genes), max(n_genes), gmt_file))
