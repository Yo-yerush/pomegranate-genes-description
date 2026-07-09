# ================================
# Download UniProt table by taxonomy ID
# taxonomy_id / organism_id: 22663
# ================================

setwd("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/pom/build_description_file_pom")

library(httr)
library(readr)
library(dplyr)

# UniProt uses organism_id for taxonomy ID queries
tax_id <- "22663"

# Requested UniProt columns
# Entry                                      = accession
# Entry Name                                 = id
# Protein names                              = protein_name
# Gene Names (primary)                       = gene_primary
# Pathway                                    = cc_pathway
# Function [CC]                              = cc_function
# EC number                                  = ec
# Gene Ontology (biological process)         = go_p
# Gene Ontology (cellular component)         = go_c
# Gene Ontology (molecular function)         = go_f
# Motif                                      = ft_motif
# Protein families                           = protein_families

fields <- c(
  "accession",
  "id",
  "protein_name",
  "gene_primary",
  "cc_pathway",
  "cc_function",
  "ec",
  "go_p",
  "go_c",
  "go_f",
  "ft_motif",
  "protein_families"
)

get_next_link <- function(link_header) {
  if (is.null(link_header)) return(NULL)
  
  m <- regmatches(
    link_header,
    regexpr("<([^>]+)>; rel=\"next\"", link_header, perl = TRUE)
  )
  
  if (length(m) == 0 || m == "") return(NULL)
  
  sub("^<([^>]+)>.*$", "\\1", m)
}

download_uniprot_paged <- function(tax_id, fields, size = 500) {
  
  base_url <- "https://rest.uniprot.org/uniprotkb/search"
  
  query_text <- paste0("organism_id:", tax_id)
  
  url <- modify_url(
    base_url,
    query = list(
      query = query_text,
      fields = paste(fields, collapse = ","),
      format = "tsv",
      size = size
    )
  )
  
  all_pages <- list()
  page <- 1
  
  while (!is.null(url)) {
    
    message("Downloading UniProt page ", page, " ...")
    
    response <- RETRY(
      verb = "GET",
      url = url,
      times = 8,
      pause_base = 2,
      pause_cap = 60,
      terminate_on = c(400, 401, 403, 404)
    )
    
    stop_for_status(response)
    
    txt <- content(response, as = "text", encoding = "UTF-8")
    
    if (grepl("Error encountered when streaming data", txt, fixed = TRUE)) {
      stop("UniProt returned a streaming error instead of data.")
    }
    
    page_df <- read_tsv(
      I(txt),
      show_col_types = FALSE,
      progress = FALSE
    )
    
    all_pages[[page]] <- page_df
    
    link_header <- headers(response)[["link"]]
    url <- get_next_link(link_header)
    
    page <- page + 1
    
    Sys.sleep(0.2)
  }
  
  bind_rows(all_pages)
}

uniprot_df <- download_uniprot_paged(
  tax_id = tax_id,
  fields = fields,
  size = 500
)

# merge to get gene_id
uni2xp <- read.csv("uni.by.XP.csv") %>%
    select(Entry, protein_id)

refseq_pom <- read.csv("pomegranate_refseq_ids.csv")

merged_ids <- merge(uni2xp, refseq_pom, by = "protein_id")

uniprot_df_out <- merge(uniprot_df, merged_ids, by = "Entry", all.y = TRUE) %>%
    arrange(protein_id) %>%
    rename(
    uniprot_id = "Entry",
    Protein_names = "Protein names",
    Symbol = "Gene Names (primary)",
    Function = "Function [CC]",
    EC_number = "EC number",
    GO_biological_process = "Gene Ontology (biological process)",
    GO_cellular_component = "Gene Ontology (cellular component)",
    GO_molecular_function = "Gene Ontology (molecular function)",
    Protein_families = "Protein families"
  ) %>%
  group_by(gene_id) %>%
  summarize(
    Symbol = gsub("LOC[0-9]+", "", first(Symbol)),
    Protein_names = first(Protein_names),
    Protein_families = first(Protein_families),
    Function = gsub("FUNCTION: ", "", first(Function)),
    Pathway = gsub("PATHWAY: ", "", first(Pathway)),
    GO_biological_process = first(GO_biological_process),
    GO_cellular_component = first(GO_cellular_component),
    GO_molecular_function = first(GO_molecular_function),
    EC_number = first(EC_number),
    Motif = first(Motif),
    protein_id = paste(unique(protein_id), collapse = "; "),
    transcript_id = paste(unique(transcript_id), collapse = "; "),
    uniprot_id = first(uniprot_id)
  )

# save csv
write.csv(
  uniprot_df_out, "pomegranate_description_file.csv", row.names = FALSE, na = "")


########################################################################################################################

## GO IDs db
GO_ids <- read.csv("C:/Users/YonatanY/Migal/Rachel Amir Team - General/yonatan/pom/RNAseq_yonatan_2021/DESeq2/uniprot/uniprot-pomegranate.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>%
    select(Entry, Gene.ontology.IDs) %>%
    rename("GO_ids" = Gene.ontology.IDs) %>%
    merge(., merged_ids, by = "Entry") %>%
    select(gene_id, GO_ids) %>%
    distinct(gene_id, .keep_all = TRUE)

write.csv(GO_ids, "pomegranate_GO_ids.csv", row.names = FALSE, na = "")












