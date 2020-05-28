# Code adapted from https://grunwaldlab.github.io/metacoder_documentation/publication--09--gene_expression.html
output_format <- "png"
output_folder <- "figures/tree_plots/"
pub_fig_folder <- "publication"
revision_n <- 1

# GO annotation file from the Gene Ontology Consortium
path_to_chicken_GO_annotation_file <- "data/chicken_go_annotation/goa_chicken.gaf"
chicken_go_annotation <- read.delim(path_to_chicken_GO_annotation_file, stringsAsFactors=FALSE, skip = 31, header = FALSE,
                                    col.names = c("DB", "DB_ID", "Symbol", "Qualifier", "GO_ID", "DB:Reference",
                                                  "Evidence_code", "With_or_from", "Aspect", "DB_Object_Name",
                                                  "DB_Object_Synonyms", "DB_Object_type", "Taxon", "Date", "Assigned_by",
                                                  "Annotation_extension", "Gene_product_form_ID"))
chicken_go_annotation <- chicken_go_annotation[chicken_go_annotation$Aspect == "P",] # Filter so only BP GOs remain

test_chicken_qlf_list <- list(qlf_chicken_UI.2,qlf_chicken_UI.4,qlf_chicken_UI.12,qlf_chicken_UI.24,
                         qlf_chicken_UI.48,qlf_chicken_UI.72)
i <- 0
while (i < length(test_chicken_qlf_list)) {
  m <- match(test_chicken_qlf_list[[i+1]]$genes$gene_name, chicken_go_annotation$Symbol)
  test_chicken_qlf_list[[i+1]]$genes$go_id <- chicken_go_annotation$GO_ID[m]
  i <- i + 1
}
test_topgenes_chicken_list <- lapply(test_chicken_qlf_list,
                                function(x) topTags(x, n = dim(test_chicken_qlf_list[[1]])[1]))

# Test variables defined
res <- test_topgenes_chicken_list[[4]]$table
output_file_name <- "chicken_tree_24h"
fdr_thresh <- 0.05
logfc_thresh <- 1
fc_scale = "both"

# Trees for all timepoints produced
test_topgenes_chicken_list_pos <- lapply(test_topgenes_chicken_list, function(x) x[x$table$logFC > 1,])
test_topgenes_chicken_list_neg <- lapply(test_topgenes_chicken_list, function(x) x[x$table$logFC < -1,])
tree_names <- c("chicken_tree_2h", "chicken_tree_4h", "chicken_tree_12h", 
                "chicken_tree_24h", "chicken_tree_48h", "chicken_tree_72h")
tree_names_pos <- c("chicken_tree_2h_pos", "chicken_tree_4h_pos", "chicken_tree_12h_pos", 
                "chicken_tree_24h_pos", "chicken_tree_48h_pos", "chicken_tree_72h_pos")
tree_names_neg <- c("chicken_tree_2h_neg", "chicken_tree_4h_neg", "chicken_tree_12h_neg", 
                "chicken_tree_24h_neg", "chicken_tree_48h_neg", "chicken_tree_72h_neg")
i <- 0
while (i < length(test_topgenes_chicken_list)) {
  i <- i + 1
  test_topgenes_chicken_list[[i]]$table$entrez_gene_id <- as.character(test_topgenes_chicken_list[[i]]$table$entrez_gene_id)
  test_topgenes_chicken_list_pos[[i]]$table$entrez_gene_id <- as.character(test_topgenes_chicken_list_pos[[i]]$table$entrez_gene_id)
  test_topgenes_chicken_list_neg[[i]]$table$entrez_gene_id <- as.character(test_topgenes_chicken_list_neg[[i]]$table$entrez_gene_id)
  create_go_tree_plot(test_topgenes_chicken_list[[i]]$table, tree_names[i], output_format = output_format, output_folder = output_folder, fc_scale = "both")
  create_go_tree_plot(test_topgenes_chicken_list_pos[[i]]$table, tree_names_pos[i], fc_scale = "pos")
  create_go_tree_plot(test_topgenes_chicken_list_neg[[i]]$table, tree_names_neg[i], fc_scale = "neg")
}

create_go_tree_plot <- function(res, output_file_name, 
                                output_format = "pdf",
                                output_folder = "results",
                                pub_fig_folder = "publication",
                                revision_n = 1,
                                fdr_thresh = 0.05,
                                logfc_thresh = 1,
                                fc_scale = "both") {
  
  result_path <- function(name) {
    file.path(output_folder, paste0(name, ".", output_format))
  }
  save_publication_fig <- function(name, figure_number) {
    file.path(result_path(name), paste0("revision_", revision_n), paste0("figure_", figure_number, "--", name, ".", output_format))
  }
  
  # Plan to improve this graph:
  #   Use the downloaded goa_chicken.gaf file to get better GO annotations for the chicken
  #   Filter the GO terms that can be mapped to by ontology, we are only interested in BP for now
  #   Increase the number of valid relationships for the GO categories in term_class
  
  library(GO.db)
  library(org.Gg.eg.db)
  
  
  # Filtering results
  
  # Remove genes with no GO annotation
  
  nrow(res)
  res <- res[!is.na(res$go_id), ]
  nrow(res)
  
  # Remove insignificant genes
  
  nrow(res)
  res <- res[(! is.na(res$FDR)) & res$FDR <= fdr_thresh, ]
  nrow(res)
  
  # Remove genes with small changes
  
  nrow(res)
  res <-  res[abs(res$logFC) >= logfc_thresh, ]
  print(nrow(res))
  a <- nrow(res) == 0
  print(a)
  if(nrow(res) == 0){
    stop("The number of rows after filtering is zero, the figure cannot be generated")
  }
  # Getting classificaiton
  
  term_class <- function(x, current = x, all_paths = TRUE, type = GOCCPARENTS, verbose = TRUE, 
                         valid_relationships = c("is_a")) {
    # Get immediate children of current taxon
    parents = tryCatch({
      possible_parents <- as.list(type[x[1]])[[1]]
      if (! is.null(valid_relationships)) {
        possible_parents <- possible_parents[names(possible_parents) %in% valid_relationships]
      }
      names(AnnotationDbi::Term(possible_parents))
    }, error = function(e) {
      c()
    })
    
    # only go down one path if desired
    if (! all_paths) {
      parents <- parents[1]
    }
    parents <- parents[parents != "all"]
    
    if (is.null(parents)) {
      return(c())
    } else if (length(parents) == 0) {
      return(paste0(collapse = "|", AnnotationDbi::Term(x)))
    } else {
      next_x <- lapply(parents, function(y) c(y, x))
      
      # Run this function on them to get their output
      child_output <- lapply(next_x, term_class, all_paths = all_paths, type = type)
      output <- unlist(child_output, recursive = FALSE)
      
      return(output)
    }
  }
  
  # cc_class <- lapply(res$go_id, term_class, all_paths = FALSE, type = GOCCPARENTS)
  # mf_class <- lapply(res$go_id, term_class, all_paths = FALSE, type = GOMFPARENTS)
  bp_class <- lapply(res$go_id, term_class, all_paths = FALSE, type = GOBPPARENTS)
  
  
  # Biological Process
  
  bp_res <- res[rep(1:nrow(res), sapply(bp_class, length)), ]
  bp_res$class <- unlist(bp_class)
  library(metacoder)
  
  obj <- parse_tax_data(as.data.frame(bp_res), class_cols = "class", class_sep = "|")
  obj$funcs <- c(obj$funcs,
                 change = function(x, subset = NULL) {
                   vapply(obs(x, "tax_data"),
                          function(i) {
                            obs_change <- obj$data$tax_data[i, ]$logFC[obj$data$tax_data[i, ]$FDR <= min_p_value]
                            mean(obs_change, na.rm = TRUE)
                          },
                          numeric(1))
                 },
                 num_changed = function(x, subset = NULL) {
                   vapply(obs(x, "tax_data"),
                          function(i) {
                            sum(obj$data$tax_data[i, ]$FDR <= min_p_value, na.rm = TRUE)
                          },
                          numeric(1))
                 })
  
  set.seed(7) #2, 4, 5, 7*, 19
  mgsub <- function(pattern, replacement, x, ...) { # from: http://stackoverflow.com/questions/15253954/replace-multiple-arguments-with-gsub
    if (length(pattern)!=length(replacement)) {
      stop("pattern and replacement do not have the same length.")
    }
    result <- x
    for (i in 1:length(pattern)) {
      result <- gsub(pattern[i], replacement[i], result, ...)
    }
    result
  }
  
  to_replace <- matrix(ncol = 2, byrow = TRUE,
                       c("regulation of growth", "",
                         "activation of innate immune response", "",
                         "system development", "",
                         "regulation of response to stimulus", "",
                         "lipid metabolic process", "",
                         "selenium compound metabolic process", "selenium metabolic process"
                       ))
  output_path <- file.path(output_folder,
                           paste0("gene_expression--biological_process",
                                  output_format))
  min_p_value <- fdr_thresh
  if (fc_scale == "pos") {
    fc_scale_interval <- c(0,10)
    fc_scale_palette <- diverging_palette()[2:3]
  } else if(fc_scale == "neg") {
    fc_scale_interval <- c(-10,0)
    fc_scale_palette <- diverging_palette()[1:2]
  } else {
    fc_scale_interval <- c(-5,5)
    fc_scale_palette <- diverging_palette()
  }
  obj %>%
    filter_taxa(num_changed > 0) %>%
    filter_taxa(n_supertaxa <= 3) %>%
    mutate_obs("plot_data", 
               taxon_id = taxon_ids,
               plotted_name = gsub("_", " ", taxon_names),
               f_change = 2^abs(change) * sign(change)) %>%
    mutate_obs("plot_data",
               short_name = vapply(plotted_name, FUN.VALUE = character(1), function(x) {
                 mgsub(pattern = to_replace[, 1], replacement =  to_replace[, 2], x, fixed = TRUE)
               })) %>%
    heat_tree(node_label = ifelse(abs(f_change) > 1, short_name, NA),
              node_size = num_changed,
              # node_size_trans = "log10",
              node_size_range = c(0.008, 0.03),
              # node_label_size_trans = "log10",
              node_label_size_range = c(0.006, 0.02),
              # edge_size_trans = "log10",
              edge_size_range = c(0.008, 0.03) / 2,
              node_color = f_change,
              node_color_trans = "linear",
              node_color_range = fc_scale_palette,
              node_color_interval = fc_scale_interval,
              edge_color_trans = "linear",
              edge_color_range = fc_scale_palette,
              edge_color_interval =  fc_scale_interval,
              node_color_axis_label = "Fold change",
              node_size_axis_label = "Number of genes",
              layout = "da", initial_layout = "re",
              output_file = result_path(output_file_name))

}











