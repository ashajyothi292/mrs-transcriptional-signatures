# Fix the gene symbols that Excel turned into dates
#
# counts_matrix_195.csv went through a spreadsheet at some point, and 27 gene
# symbols came out as dates (MARCH1 -> "1-Mar", SEPT2 -> "2-Sep", DEC1 -> "1-Dec").
# Two of them even ended up duplicated, because MARCH1 and MARC1 both become
# "1-Mar" (same for MARCH2 / MARC2), which stops read.csv(row.names = 1) working.
# This is the problem described in Ziemann et al. 2016 (Genome Biology).
#
# The rows are in genome order, so the duplicates can be told apart by their
# neighbours: the first "2-Mar"/"1-Mar" pair sits on chromosome 1 next to
# C1orf115 (that is MARC2 and MARC1), the later ones sit next to TKTL2 on
# chromosome 4 (MARCH1) and RAB11B on chromosome 19 (MARCH2).
#
# I restore them to the current HGNC names (MARCHF1, SEPTIN2, MTARC1 ...) because
# that is what the MSigDB gene sets use.
#
# Run from the repository root:  Rscript scripts/02_fix_gene_symbols.R

counts <- read.csv("data/raw/counts_matrix_195.csv", check.names = FALSE)
genes <- counts[[1]]

date_rows <- grep("^[0-9]{1,2}-(Mar|Sep|Dec)$", genes)
cat("Date-converted symbols found:", length(date_rows), "\n")

fixed <- genes
seen_mar <- c()
for (i in date_rows) {
  num <- sub("-.*", "", genes[i])
  mon <- sub(".*-", "", genes[i])
  if (mon == "Sep") {
    fixed[i] <- paste0("SEPTIN", num)
  } else if (mon == "Dec") {
    fixed[i] <- paste0("DELEC", num)          # DEC1 is now called DELEC1
  } else if (num %in% c("1", "2") && !(num %in% seen_mar)) {
    fixed[i] <- paste0("MTARC", num)          # first 1-Mar / 2-Mar = MARC1 / MARC2 (chr1)
    seen_mar <- c(seen_mar, num)
  } else {
    fixed[i] <- paste0("MARCHF", num)
  }
}
print(data.frame(row = date_rows, before = genes[date_rows], after = fixed[date_rows]))

stopifnot(!any(duplicated(fixed)))
counts[[1]] <- fixed
names(counts)[1] <- ""
write.csv(counts, "data/raw/counts_matrix_195_fixed.csv", row.names = FALSE)
cat(sprintf("\nWrote data/raw/counts_matrix_195_fixed.csv: %d genes x %d samples\n",
            nrow(counts), ncol(counts) - 1))
