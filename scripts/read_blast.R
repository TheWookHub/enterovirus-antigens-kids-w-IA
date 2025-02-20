# Function to read BLAST output
read_blast <- function(input_blast_result) {
  read_tsv(input_blast_result, col_names = c("qaccver", "saccver", "pident", "nident", "length", "evalue", "bitscore", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "qseq", "sseq", "ppos", "stitle", "frames"))
}