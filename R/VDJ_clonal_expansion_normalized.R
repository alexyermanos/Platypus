# the possible normalization method is: split_capture
# the variable groups_of_capture is a list of lists indicating which mice (so the FB assignments) belong to which capture

#--------------------------------------
## Example code

# # Define capture groups (you can add here also the not assigned cells)
# capture1 <- list("Mouse 1","Mouse 2", "Mouse 3", "Not assignable_1")
# capture2 <- list( "Mouse 4", "Not assignable_2")
# 
# # for the until_which_clone_normalize we want the number of clones which have minimum 2 cells. this will be the variable: clones_with_two_cells
# clono_df <- as.data.frame(table(VGM$clonotype_id))
# clono_df <- clono_df[order(clono_df$Freq, 
#                            decreasing = T), ]
# clono_df <- clono_df[which(clono_df$Freq > 1),]
# clones_with_two_cells <- nrow(clono_df)
# 
# clonal_out <- VDJ_clonal_expansion_normalized(VDJ = VGM, clones = "30", group.by = "sample_id", color.by = "isotype", isotypes.to.plot = "all", species = "Mouse", treat.incomplete.clones = "exclude", treat.incomplete.cells = "proportional", normalization_method = "split_capture", groups_of_capture = list(capture1, capture2), until_which_clone_normalize = clones_with_two_cells)
# 

#--------------------------------------


VDJ_clonal_expansion_normalized <- function (VDJ, celltype, clones, subtypes, isotypes.to.plot, 
                                                               species, treat.incomplete.clones, treat.incomplete.cells, 
                                                               group.by, color.by, variant.plot, normalization_method, groups_of_capture, list_of_mice, until_which_clone_normalize) 
{
  
  
  VDJ.matrix <- VDJ
  VDJ <- NULL
  ClonalRank <- NULL
  Counts <- NULL
  sum_counts <- NULL
  Isotype <- NULL
  Color <- NULL
  pasted_variants <- NULL
  isotype <- NULL
  colors <- NULL
  match_ex_crit <- NULL
  clonotype_id <- NULL
  variant <- NULL
  n <- NULL
  if (missing(clones)) 
    clones <- 50
  if (missing(subtypes)) 
    subtypes <- FALSE
  if (missing(isotypes.to.plot)) 
    isotypes.to.plot <- "all"
  if (missing(species)) 
    species <- "Human"
  if (missing(treat.incomplete.cells)) 
    treat.incomplete.cells <- "proportional"
  if (missing(treat.incomplete.clones)) 
    treat.incomplete.clones <- "exclude"
  if (!treat.incomplete.cells %in% c("exclude", "proportional")) {
    stop("Please set treat.incomplete.cells to either 'proportional' or 'exclude'. 'Proportional' will assign cells of a clonotype missing a VDJ chain proportionally to the isotypes present in that clone (Default). 'exclude' will remove all cells missing a VDJ chain and thereby also alter clonotype frequencies")
  }
  if (!treat.incomplete.clones %in% c("exclude", "include")) {
    stop("Please set treat.incomplete.clones to either 'include' or 'exclude'. 'include' will show also clones which of which no cell has a VDJ chain. 'exclude' will remove such clones (default)")
  }
  platypus.version <- "v3"
  if (missing(celltype)) {
    if (stringr::str_detect(VDJ.matrix$celltype[1], "B")) {
      celltype = "Bcells"
    }
    else if (stringr::str_detect(VDJ.matrix$celltype[1], 
                                 "T")) {
      celltype = "Tcells"
    }
    else {
      stop("No celltype found in celltype column. celltype column must contain either 'B cells' or 'T cells'")
    }
  }
  if (missing(group.by)) {
    group.by <- "sample_id"
  }
  if (missing(variant.plot)) {
    variant.plot <- FALSE
  }
  if (group.by == "none") {
    group.by <- "ungroup"
    VDJ.matrix$ungroup <- 1
  }
  else if (!group.by %in% names(VDJ.matrix) & group.by != 
           "sample_id" & group.by == "none") {
    stop("Please provide a valid column name of the VDJ.matrix in group.by")
  }
  if (missing(color.by)) {
    color.by <- "isotype"
  }
  if (missing(normalization_method)) {
    normalization_method <- NULL
  }
  if (missing(groups_of_capture)) {
    groups_of_capture <- NULL
  }
  if (missing(list_of_mice)) {
    list_of_mice <- NULL
  }
  if (missing(until_which_clone_normalize)) {
    until_which_clone_normalize <- 100
  }
  if (!color.by %in% names(VDJ.matrix) & color.by != "isotype") {
    stop("Please provide a valid column name of the VDJ.matrix in color.by, or set to 'isotype' for default (also for Tcells)")
  }
  else if (color.by %in% names(VDJ.matrix) & color.by != "isotype") {
    unique_colors <- grDevices::rainbow(length(unique(VDJ.matrix[, 
                                                                 color.by])))
    names(unique_colors) <- unique(VDJ.matrix[, color.by])
  }
  clones <- strtoi(clones)
  
  #for capture normalization
  if (normalization_method == "split_capture" & is.null(groups_of_capture) == FALSE){
    # get the values of the numbers of cells in each capture
    counts_of_each_mouse <- VDJ.matrix %>% count(FB_assignment)
    number_of_mice <- length(counts_of_each_mouse$n)
    # figure out how many mice in each capture
    n_capture1 <- length(groups_of_capture[[1]])
    n_capture2 <- length(groups_of_capture[[2]])
    # loop through capture 1
    total_n_capture1 <- 0
    for (m in 1:n_capture1){
      # OLD CODE FROM ME
      #name_of_mouse <- groups_of_capture[[1]][[m]]
      # get the number of this mouse and add it to the variable
      #total_n_capture1 <- total_n_capture1 + counts_of_each_mouse[which(counts_of_each_mouse$FB_assignment == name_of_mouse),]$n
      # NEW CODE FROM ME
      name_of_mouse <- groups_of_capture[[1]][[m]]
      filtered_for_mouse <- subset(VDJ.matrix, FB_assignment == name_of_mouse)
      sum <- 0
      unique_clonotypes <- unique(filtered_for_mouse$clonotype_id)
      for (specific_clonotype in unique_clonotypes) {
        value_frequency <- unique(subset(filtered_for_mouse, clonotype_id == specific_clonotype)$clonotype_frequency)
        sum <- sum + value_frequency
      }
      total_n_capture1 <- total_n_capture1 + sum
    }
    # loop through capture 2
    total_n_capture2 <- 0
    for (m in 1:n_capture2){
      if (n_capture2 != 0) {
        # OLD CODE FROM ME
        #name_of_mouse <- groups_of_capture[[2]][[m]]
        ## get the number of this mouse and add it to the variable
        #total_n_capture2 <- total_n_capture2 + counts_of_each_mouse[which(counts_of_each_mouse$FB_assignment == name_of_mouse),]$n
        # NEW CODE FROM ME
        name_of_mouse <- groups_of_capture[[2]][[m]]
        filtered_for_mouse <- subset(VDJ.matrix, FB_assignment == name_of_mouse)
        sum <- 0
        unique_clonotypes <- unique(filtered_for_mouse$clonotype_id)
        for (specific_clonotype in unique_clonotypes) {
          value_frequency <- unique(subset(filtered_for_mouse, clonotype_id == specific_clonotype)$clonotype_frequency)
          sum <- sum + value_frequency
        }
        total_n_capture2 <- total_n_capture2 + sum
      }
    }
    print(total_n_capture1)
    print(total_n_capture2)
  }
  
  #for mouse normalization
  if (normalization_method == "split_mouse" & is.null(list_of_mice) == FALSE){
    # get the values of the numbers of cells in each mouse
    counts_of_each_mouse <- VDJ.matrix %>% count(FB_assignment)
    # how many mice do you have
    number_of_mice <- length(list_of_mice)
  }
  
  
  VDJ.matrix$clonotype_frequency <- VDJ.matrix$clonotype_frequency/1000
  if (celltype == "Bcells") {
    if (!inherits(VDJ.matrix, "data.frame")) {
      stop("Please provide a VDJ dataframe. (VDJ_GEX_matrix_out[[1]]")
    }
    VDJ.matrix[, group.by] <- as.character(VDJ.matrix[, 
                                                      group.by])
    sample.names <- unique(VDJ.matrix[, group.by])
    sample.names[is.na(sample.names)] <- "NONE"
    VDJ.matrix[is.na(VDJ.matrix[, group.by]), group.by] <- "NONE"
    VDJ.matrix.all <- VDJ.matrix
    VDJ_per_clone_output_all <- list()
    clones_per_isotype <- list()
    clones_per_isotype_all <- list()
    output_plot <- list()
    variant_df_list <- list()
    if (subtypes == FALSE | color.by != "isotype") {
      if (variant.plot == FALSE) {
        for (i in 1:length(sample.names)) {
          VDJ.matrix <- subset(VDJ.matrix.all, VDJ.matrix.all[, 
                                                              group.by] == sample.names[i])
          clono_freq <- as.data.frame(table(VDJ.matrix$clonotype_id))
          clono_freq <- clono_freq[order(clono_freq$Freq, 
                                         decreasing = T), ]
          if (is.null(VDJ.matrix$overlap_type)) {
            curr_rep_iso <- VDJ.matrix[, c("barcode", 
                                           "clonotype_id", "VDJ_cgene", "VDJ_cdr3s_aa", 
                                           "VJ_cdr3s_aa", "FB_assignment")]
          } else {
            curr_rep_iso <- VDJ.matrix[, c("barcode", 
                                           "clonotype_id", "VDJ_cgene", "VDJ_cdr3s_aa", 
                                           "VJ_cdr3s_aa", "FB_assignment", "overlap_type")]
          }
          names(curr_rep_iso)[3] <- "isotype"
          if (color.by != "isotype") {
            curr_rep_iso$colors <- as.character(VDJ.matrix[, 
                                                           color.by])
            curr_rep_iso$colors[which(is.na(curr_rep_iso$colors))] <- "None"
            if (inherits(VDJ.matrix[, color.by], "factor")) {
              curr_rep_iso$colors <- ordered(as.factor(curr_rep_iso$colors), 
                                             levels = c(levels(VDJ.matrix[, color.by]), 
                                                        "None"))
            }
            else {
              curr_rep_iso$colors <- ordered(as.factor(curr_rep_iso$colors), 
                                             levels = c(as.character(unique(VDJ.matrix[, 
                                                                                       color.by])), "None"))
            }
          }
          if (treat.incomplete.clones == "exclude") {
            clones_to_del <- c()
            for (k in 1:nrow(clono_freq)) {
              if (all(curr_rep_iso[which(curr_rep_iso$clonotype_id == 
                                         clono_freq[k, 1]), "isotype"] == "")) {
                clones_to_del <- append(clones_to_del, 
                                        k)
              }
            }
            if (length(clones_to_del > 0)) {
              clono_freq <- clono_freq[-clones_to_del, 
              ]
            }
            clono_freq <- clono_freq[order(clono_freq$Freq, 
                                           decreasing = T), ]
          }
          clones_per_isotype <- list()
          # just take in consideration the first 500 clones
          for (j in 1:until_which_clone_normalize) {
            curr_clone <- curr_rep_iso[which(curr_rep_iso$clonotype_id == 
                                               clono_freq[j, 1]), ]
            if (nrow(curr_clone) > 0) {
              if (color.by == "isotype") {
                curr_clone$isotype <- stringr::str_split(curr_clone$isotype, 
                                                         ";", simplify = T)[, 1]
                if (treat.incomplete.cells == "proportional" & 
                    stringr::str_detect(paste0(curr_clone$isotype, 
                                               collapse = ";"), pattern = "IG") == 
                    T) {
                  props <- table(curr_clone$isotype[which(stringr::str_detect(curr_clone$isotype, 
                                                                              pattern = "IG") == T)])
                  n_total <- nrow(curr_clone)
                  props <- round(props/sum(props) * 
                                   n_total, 0)
                  if (sum(props) > n_total) {
                    props[which.max(props)] <- props[which.max(props)] - 
                      1
                  }
                  else if (sum(props) < n_total) {
                    props[which.min(props)] <- props[which.min(props)] + 
                      1
                  }
                  curr_clone$isotype <- rep.int(names(props), 
                                                props)
                }
                else if (treat.incomplete.cells == "proportional" & 
                         stringr::str_detect(paste0(curr_clone$isotype, 
                                                    collapse = ";"), pattern = "IG") == 
                         F) {
                  curr_clone$isotype <- "None"
                }
                clones_per_isotype[[j]] <- data.frame(Counts = rep(0, 
                                                                   6), Color = rep("", 6), Isotype = rep("", 
                                                                                                         6), ClonalRank = rep("", 6), clonotype_id = rep(curr_clone$clonotype_id[1], 
                                                                                                                                                         6), VDJ_cdr3s_aa = rep(curr_clone$VDJ_cdr3s_aa[which(curr_clone$VDJ_cdr3s_aa != 
                                                                                                                                                                                                                "")][1], 6), VJ_cdr3s_aa = rep(curr_clone$VJ_cdr3s_aa[which(curr_clone$VJ_cdr3s_aa != 
                                                                                                                                                                                                                                                                              "")][1], 6), FB_assignment = rep("",6), barcode = rep(paste0(curr_clone$barcode, 
                                                                                                                                                                                                                                                                                                                                           collapse = ";"), 6))
                
                
                
                # from me 
                if (normalization_method == "split_capture") {
                  
                  clones_per_isotype[[j]]$FB_assignment <- curr_clone$FB_assignment[1]
                  
                  if (clones_per_isotype[[j]]$FB_assignment[1] %in% groups_of_capture[[1]]) {
                    
                    clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHG"))/total_n_capture1
                    clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHM"))/total_n_capture1
                    clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHA"))/total_n_capture1
                    clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHD"))/total_n_capture1
                    clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHE"))/total_n_capture1
                    clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "None"))/total_n_capture1
                  }
                  else if (clones_per_isotype[[j]]$FB_assignment[1] %in% groups_of_capture[[2]]) {
                    
                    clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHG"))/total_n_capture2
                    clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHM"))/total_n_capture2
                    clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHA"))/total_n_capture2
                    clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHD"))/total_n_capture2
                    clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHE"))/total_n_capture2
                    clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "None"))/total_n_capture2
                  }
                  else if (clones_per_isotype[[j]]$FB_assignment[1] == "Not assignable") {
                    
                    
                    
                    # for the case that the capture 2 is not empty, just divide by mean
                    
                    clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHG"))/((total_n_capture1+total_n_capture2)/2)
                    clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHM"))/((total_n_capture1+total_n_capture2)/2)
                    clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHA"))/((total_n_capture1+total_n_capture2)/2)
                    clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHD"))/((total_n_capture1+total_n_capture2)/2)
                    clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHE"))/((total_n_capture1+total_n_capture2)/2)
                    clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "None"))/((total_n_capture1+total_n_capture2)/2)
                    
                  }
                  
                  
                  
                  
                  
                }
                
                
                
                
                
                # from me 
                else if (normalization_method == "split_mouse") {
                  if (length(unique(curr_clone$FB_assignment)) == 1) {
                    clones_per_isotype[[j]]$FB_assignment <- curr_clone$FB_assignment[1]
                    
                    for (mouse in list_of_mice) {
                      if (clones_per_isotype[[j]]$FB_assignment[1] == mouse) {
                        counts_of_this_mouse <- counts_of_each_mouse[which(counts_of_each_mouse$FB_assignment == mouse),]$n
                        
                        clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                    "IGHG"))/counts_of_this_mouse
                        clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                    "IGHM"))/counts_of_this_mouse
                        clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                    "IGHA"))/counts_of_this_mouse
                        clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                    "IGHD"))/counts_of_this_mouse
                        clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                    "IGHE"))/counts_of_this_mouse
                        clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                    "None"))/counts_of_this_mouse
                      }
                      
                    }
                  }
                  else if (length(unique(curr_clone$FB_assignment)) > 1) {
                    mice_in_this_clone <- unique(curr_clone$FB_assignment)
                    factor_quotient <- 0
                    total_counts_all_mice <- sum(table(curr_clone$FB_assignment))
                    for (mouse_in_clone in mice_in_this_clone){
                      occurences <- table(curr_clone$FB_assignment)[mouse_in_clone][[1]]
                      counts_mouse <- counts_of_each_mouse[which(counts_of_each_mouse$FB_assignment == mouse_in_clone),]$n
                      factor_quotient <- factor_quotient + occurences/counts_mouse
                    }
                    to_divide_by <- total_counts_all_mice/factor_quotient
                    
                    clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHG"))/to_divide_by
                    clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHM"))/to_divide_by
                    clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHA"))/to_divide_by
                    clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHD"))/to_divide_by
                    clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "IGHE"))/to_divide_by
                    clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                                "None"))/to_divide_by
                  }
                }
                
                
                
                if (color.by == "isotype") {
                  clones_per_isotype[[j]]$Color <- c("green4", 
                                                     "black", "red", "blue", "purple", 
                                                     "gray")
                  clones_per_isotype[[j]]$Isotype <- c("IGHG", 
                                                       "IGHM", "IGHA", "IGHD", "IGHE", 
                                                       "Unknown")
                }
                else {
                  clones_per_isotype[[j]]$Isotype <- names(which.max(table(curr_clone$colors)))
                }
                clones_per_isotype[[j]]$ClonalRank <- j
                
              }
              else if (color.by != "isotype") {
                color_cur_clone <- unique(curr_clone$colors)
                n_color_cur_clone <- length(unique(curr_clone$colors))
                clones_per_isotype[[j]] <- data.frame(Counts = rep(0, 
                                                                   n_color_cur_clone), Color = color_cur_clone, 
                                                      #from me
                                                      FB_assignment = rep(curr_clone$FB_assignment[1], n_color_cur_clone),
                                                      ClonalRank = rep("", n_color_cur_clone), 
                                                      clonotype_id = rep(curr_clone$clonotype_id[1], 
                                                                         n_color_cur_clone), VDJ_cdr3s_aa = rep(curr_clone$VDJ_cdr3s_aa[which(curr_clone$VDJ_cdr3s_aa != 
                                                                                                                                                "")][1], n_color_cur_clone), VJ_cdr3s_aa = rep(curr_clone$VJ_cdr3s_aa[which(curr_clone$VJ_cdr3s_aa != 
                                                                                                                                                                                                                              "")][1], n_color_cur_clone), barcode = rep(paste0(curr_clone$barcode, 
                                                                                                                                                                                                                                                                                collapse = ";"), n_color_cur_clone))
                #from me
                #clones_per_isotype[[j]]$FB_assignment <- curr_clone$FB_assignment[1]
                for (k in 1:nrow(clones_per_isotype[[j]])) {
                  
                  if (color.by == 'FB_assignment') {
                    
                    
                    if (normalization_method == 'split_capture') {
                      
                      
                      if (clones_per_isotype[[j]]$Color[k] %in% groups_of_capture[[1]]) {
                        clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                                       paste0(curr_clone$colors, collapse = "/ /"), 
                                                                                       "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                              "/"))/total_n_capture1
                      }
                      else if (clones_per_isotype[[j]]$Color[k] %in% groups_of_capture[[2]]) {
                        clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                                       paste0(curr_clone$colors, collapse = "/ /"), 
                                                                                       "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                              "/"))/total_n_capture2
                      }
                      else if (clones_per_isotype[[j]]$Color[k] == "Not assignable") {
                        
                        clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                                       paste0(curr_clone$colors, collapse = "/ /"), 
                                                                                       "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                              "/"))/((total_n_capture1+total_n_capture2)/2)
                      }
                      
                    }
                    else if (normalization_method == 'split_mouse') {
                      
                      for (mouse_FB in list_of_mice) {
                        if (clones_per_isotype[[j]]$Color[k] == mouse_FB) {
                          counts_of_this_mouse <- counts_of_each_mouse[which(counts_of_each_mouse$FB_assignment == mouse_FB),]$n
                          
                          clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                                         paste0(curr_clone$colors, collapse = "/ /"), 
                                                                                         "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                                "/"))/counts_of_this_mouse
                        }
                      }
                    }
                    
                    
                  }
                  
                  
                  
                  else if (color.by == 'overlap_type') {
                    
                    if (normalization_method == 'split_capture') {
                      
                      
                      if (clones_per_isotype[[j]]$FB_assignment[k] %in% groups_of_capture[[1]]) {
                        clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                                       paste0(curr_clone$colors, collapse = "/ /"), 
                                                                                       "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                              "/"))/total_n_capture1
                      }
                      else if (clones_per_isotype[[j]]$FB_assignment[k] %in% groups_of_capture[[2]]) {
                        clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                                       paste0(curr_clone$colors, collapse = "/ /"), 
                                                                                       "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                              "/"))/total_n_capture2
                      }
                      else if (clones_per_isotype[[j]]$FB_assignment[k] == "Not assignable") {
                        
                        clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                                       paste0(curr_clone$colors, collapse = "/ /"), 
                                                                                       "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                              "/"))/((total_n_capture1+total_n_capture2)/2)
                      }
                    }
                    
                  }
                  
                }
                clones_per_isotype[[j]]$ClonalRank <- j
              }
            }
          }
          clones_per_isotype_all[[i]] <- do.call("rbind", 
                                                 clones_per_isotype)
          if (treat.incomplete.cells == "exclude" & 
              color.by == "isotype") {
            rank_raw <- as.data.frame(clones_per_isotype_all[[i]] %>% 
                                        dplyr::group_by(ClonalRank) %>% dplyr::summarise(sum_counts = sum(Counts)) %>% 
                                        dplyr::arrange(dplyr::desc(sum_counts)) %>% 
                                        dplyr::mutate(rank = 1:length(unique(ClonalRank))))
            clones_per_isotype_all[[i]]$ClonalRank_2 <- 0
            for (l in 1:nrow(rank_raw)) {
              clones_per_isotype_all[[i]]$ClonalRank_2[which(clones_per_isotype_all[[i]]$ClonalRank == 
                                                               rank_raw$ClonalRank[l])] <- rank_raw$rank[l]
            }
            clones_per_isotype_all[[i]]$ClonalRank <- clones_per_isotype_all[[i]]$ClonalRank_2
            message("New ranking based only on present VDJ chains: ")
            message(unique(clones_per_isotype_all[[i]]$clonotype_id))
          }
          if (color.by == "isotype" & isotypes.to.plot[1] != 
              "all") {
            if (!any(isotypes.to.plot %in% unique(clones_per_isotype_all[[i]]$Isotype))) {
              stop("isotype.to.plot input not found in dataframe. Please check if the isotype is spelled correctly")
            }
            to_del <- c()
            for (k in 1:length(unique(clones_per_isotype_all[[i]]$ClonalRank))) {
              if (sum(clones_per_isotype_all[[i]]$Counts[which(clones_per_isotype_all[[i]]$ClonalRank == 
                                                               unique(clones_per_isotype_all[[i]]$ClonalRank)[k] & 
                                                               clones_per_isotype_all[[i]]$Isotype %in% 
                                                               isotypes.to.plot)]) == 0) {
                to_del <- append(to_del, unique(clones_per_isotype_all[[i]]$ClonalRank)[k])
              }
            }
            clones_per_isotype_all[[i]] <- subset(clones_per_isotype_all[[i]], 
                                                  !ClonalRank %in% to_del)
            rank_raw <- as.data.frame(clones_per_isotype_all[[i]] %>% 
                                        dplyr::group_by(ClonalRank) %>% dplyr::summarise(sum_counts = sum(Counts)) %>% 
                                        dplyr::arrange(dplyr::desc(sum_counts)) %>% 
                                        dplyr::mutate(rank = 1:length(unique(ClonalRank))))
            clones_per_isotype_all[[i]]$ClonalRank_2 <- 0
            for (l in 1:nrow(rank_raw)) {
              clones_per_isotype_all[[i]]$ClonalRank_2[which(clones_per_isotype_all[[i]]$ClonalRank == 
                                                               rank_raw$ClonalRank[l])] <- rank_raw$rank[l]
            }
            clones_per_isotype_all[[i]]$ClonalRank <- clones_per_isotype_all[[i]]$ClonalRank_2
            message("New ranking based only on selected isotypes: ")
            message(unique(clones_per_isotype_all[[i]]$clonotype_id))
          }
          if (color.by == "isotype") {
            
            
            
            
            
            
            #loop through normalization methods
            # if (normalization_method == "split_capture") {
            #   
            #   # loop through every row of the clones_per_isotype object and see which mouse the clone belongs to
            #   for (row in 1:nrow(clones_per_isotype_all[[i]])) {
            #     #check if mouse is in capture 1
            #     if (clones_per_isotype_all[[i]]$FB_assignment[row] %in% groups_of_capture[[1]]){
            #       clones_per_isotype_all[[i]]$Counts[row] <- clones_per_isotype_all[[i]]$Counts[row]/total_n_capture1
            #     }
            #     #check if mouse is in capture 2
            #     else if (clones_per_isotype_all[[i]]$FB_assignment[row] %in% groups_of_capture[[2]]){
            #       clones_per_isotype_all[[i]]$Counts[row] <- clones_per_isotype_all[[i]]$Counts[row]/total_n_capture2
            #     }
            #     # if clone is not assignable
            #     else if (clones_per_isotype_all[[i]]$FB_assignment[row] == "Not assignable"){
            #       clones_per_isotype_all[[i]]$Counts[row] <- clones_per_isotype_all[[i]]$Counts[row]/total_n_capture1
            #     }
            #   }
            #   # sort the data frame according to the new Counts values
            #   clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][with(clones_per_isotype_all[[i]], order((-clones_per_isotype_all[[i]]$Counts))),]
            #   # just take the top 30, 50, or 100 clones(according to what is the "clones" value)
            #   clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][1:clones, ]
            #   
            #   # we have to delete rows of clones that have e.g. one part IgM and one part IgA (because they add up and falsify the plot)
            #   # loop through the data frame
            #   for (row1 in 1:nrow(clones_per_isotype_all[[i]])) {
            #     # loop again
            #     for (row2 in 1:nrow(clones_per_isotype_all[[i]])) {
            #       # check if clonalranks are the same
            #       if ((is.na(clones_per_isotype_all[[i]]$ClonalRank[row1]) == FALSE) & (is.na(clones_per_isotype_all[[i]]$ClonalRank[row2]) == FALSE) & (clones_per_isotype_all[[i]]$ClonalRank[row1] == clones_per_isotype_all[[i]]$ClonalRank[row2]) & (row1 != row2)) {
            #         # check which one has more expansion
            #         if (clones_per_isotype_all[[i]]$Counts[row1] > clones_per_isotype_all[[i]]$Counts[row2]) {
            #           clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][-c(row2),] }
            #           else {
            #             clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][-c(row1),]
            #           }
            #         
            #       }
            #     }
            #   }
            # }  
            
            # sort the data frame according to the new Counts values
            sorted_df <- clones_per_isotype_all[[i]][with(clones_per_isotype_all[[i]], order((-clones_per_isotype_all[[i]]$Counts))),]
            # get the Clonal Ranks with the highest Counts
            #just take the top 30, 50, or 100 clones(according to what is the "clones" value)
            final_clonal_ranks <- sorted_df$ClonalRank[1:clones]
            clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][which(clones_per_isotype_all[[i]]$ClonalRank %in% final_clonal_ranks),]    
            
            
            
            # this one is changed
            output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], 
                                                ggplot2::aes(fill = Isotype, 
                                                             x = reorder(ClonalRank, -Counts), y = Counts)) + ggplot2::geom_bar(stat = "identity", 
                                                                                                                                width = 0.6, color = "black") + ggplot2::theme_bw() + 
              ggplot2::scale_fill_manual(values = c(IGHG = "green4", 
                                                    IGHM = "black", IGHA = "red3", IGHD = "blue", 
                                                    IGHE = "purple", Unknown = "gray")) + 
              ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + 
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::scale_x_discrete(expand = c(0, 
                                                                                                0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                      x = "Clonal rank", y = "Number of cells")
            
            
            
          }
          else {
            
            
            
            
            # sort the data frame according to the new Counts values
            sorted_df <- clones_per_isotype_all[[i]][with(clones_per_isotype_all[[i]], order((-clones_per_isotype_all[[i]]$Counts))),]
            # get the Clonal Ranks with the highest Counts
            #just take the top 30, 50, or 100 clones(according to what is the "clones" value)
            
            # CHANGE THIS following DO A LOOP WHERE YOU ONLY ADD CLONES TO FINAL CLONAL RANKS IF THEY ARE NOT IN THIS LIST YET AND THE END LIST SHOULD BE LENGTH 30, SO YOU CAN DO WHILE LENGTH LESS THAN 30
            final_clonal_ranks <- vector(mode='list', length=0)
            num_clone = 1
            numbering = 1
            while (length(final_clonal_ranks) < clones) {
              if (sorted_df$ClonalRank[num_clone] %in% final_clonal_ranks == FALSE){
                final_clonal_ranks[numbering] <- sorted_df$ClonalRank[num_clone]
                #final_clonal_ranks <- append(final_clonal_ranks, sorted_df$ClonalRank[num_clone])
                numbering = numbering + 1
              }
              num_clone = num_clone + 1
            }
            clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][which(clones_per_isotype_all[[i]]$ClonalRank %in% final_clonal_ranks),]    
            
            # now a new column in the dataframe has to be made combined_counts, which adds all the counts of one clone. this is needed for the plot afterwards, to sort correctly
            # loop through list with clones
            clones_per_isotype_all[[i]][ , 'combined_counts'] <- NA
            for (single_clone in final_clonal_ranks){
              for (row_df in 1:nrow(clones_per_isotype_all[[i]])) {
                if (clones_per_isotype_all[[i]]$ClonalRank[row_df] == single_clone) {
                  # create a dataframe which only contains rows with this clone
                  filtered_df <- clones_per_isotype_all[[i]][which(clones_per_isotype_all[[i]]$ClonalRank == single_clone),]    
                  # add all the counts together of this clone
                  clones_per_isotype_all[[i]]$combined_counts[row_df] <- sum(filtered_df$Counts)
                }
              }
            }
            
            #loop through normalization methods
            # if (normalization_method == "split_capture") {
            #   
            #   # loop through every row of the clones_per_isotype object and see which mouse the clone belongs to
            #   for (row in 1:nrow(clones_per_isotype_all[[i]])) {
            #     #check if mouse is in capture 1
            #     if (clones_per_isotype_all[[i]]$FB_assignment[row] %in% groups_of_capture[[1]]){
            #       clones_per_isotype_all[[i]]$Counts[row] <- clones_per_isotype_all[[i]]$Counts[row]/total_n_capture1
            #     }
            #     #check if mouse is in capture 2
            #     else if (clones_per_isotype_all[[i]]$FB_assignment[row] %in% groups_of_capture[[2]]){
            #       clones_per_isotype_all[[i]]$Counts[row] <- clones_per_isotype_all[[i]]$Counts[row]/total_n_capture2
            #     }
            #     # if clone is not assignable
            #     else if (clones_per_isotype_all[[i]]$FB_assignment[row] == "Not assignable"){
            #       clones_per_isotype_all[[i]]$Counts[row] <- clones_per_isotype_all[[i]]$Counts[row]/total_n_capture1
            #     }
            #   }
            #   # sort the data frame according to the new Counts values
            #   clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][with(clones_per_isotype_all[[i]], order((-clones_per_isotype_all[[i]]$Counts))),]
            #   # just take the top 30, 50, or 100 clones(according to what is the "clones" value)
            #   clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][1:clones, ]
            #   
            #   # we have to delete rows of clones that have e.g. one part IgM and one part IgA (because they add up and falsify the plot)
            #   # loop through the data frame
            #   for (row1 in 1:nrow(clones_per_isotype_all[[i]])) {
            #     # loop again
            #     for (row2 in 1:nrow(clones_per_isotype_all[[i]])) {
            #       # check if clonalranks are the same
            #       if ((is.na(clones_per_isotype_all[[i]]$ClonalRank[row1]) == FALSE) & (is.na(clones_per_isotype_all[[i]]$ClonalRank[row2]) == FALSE) & (clones_per_isotype_all[[i]]$ClonalRank[row1] == clones_per_isotype_all[[i]]$ClonalRank[row2]) & (row1 != row2)) {
            #         # check which one has more expansion
            #         if (clones_per_isotype_all[[i]]$Counts[row1] > clones_per_isotype_all[[i]]$Counts[row2]) {
            #           clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][-c(row2),] }
            #           else {
            #             clones_per_isotype_all[[i]] <- clones_per_isotype_all[[i]][-c(row1),]
            #           }
            #         
            #       }
            #     }
            #   }
            # }  
            
            
            
            
            # this one is changed too
            output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], 
                                                ggplot2::aes(fill = Color, x = reorder(ClonalRank, -combined_counts), y = Counts)) + ggplot2::geom_bar(stat = "identity", 
                                                                                                                                                       width = 0.6, color = "black") + ggplot2::theme_bw() + 
              ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + 
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::scale_x_discrete(expand = c(0, 
                                                                                                0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                      x = "Clonal rank", y = "Number of cells", 
                                                                                                                      fill = color.by) + ggplot2::scale_fill_manual(values = grDevices::rainbow(n = length(unique(clones_per_isotype_all[[i]]$Color))))
          }
        }
        names(clones_per_isotype_all) <- sample.names
        return(list(output_plot, clones_per_isotype_all))
      }
      if (variant.plot == TRUE) {
        for (i in 1:length(sample.names)) {
          message(paste0("Starting sample ", i, "/", 
                         length(sample.names)))
          VDJ.matrix <- subset(VDJ.matrix.all, VDJ.matrix.all[, 
                                                              group.by] == sample.names[i])
          s_id = VDJ.matrix$clonotype_id
          s_id_num = readr::parse_number(s_id)
          s_id_num = sprintf("%05d", s_id_num)
          s_id_num = paste0("clonotype", s_id_num)
          VDJ.matrix$clonotype_id <- s_id_num
          s_id = NULL
          s_id_num = NULL
          VDJ.matrix$clonotype_id = paste0("old", VDJ.matrix$clonotype_id)
          rank_raw <- as.data.frame(VDJ.matrix %>% dplyr::group_by(clonotype_id) %>% 
                                      dplyr::summarise(n = n()))
          rank_raw <- rank_raw[order(rank_raw$n, rank_raw$clonotype_id, 
                                     decreasing = T), ]
          rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
          rank_raw$new_clonotype = sprintf("%05d", rank_raw$new_clonotype)
          rank_raw$new_clonotype = paste0("new", rank_raw$new_clonotype)
          for (z in 1:nrow(rank_raw)) {
            VDJ.matrix$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                            rank_raw$new_clonotype[z], VDJ.matrix$clonotype_id)
          }
          VDJ.matrix$clonotype_id <- paste0("clonotype", 
                                            sprintf("%05d", readr::parse_number(VDJ.matrix$clonotype_id)))
          clono_freq <- as.data.frame(table(VDJ.matrix$clonotype_id))
          clono_freq <- clono_freq[order(clono_freq$Freq, 
                                         decreasing = T), ]
          is.trimmed = VDJ.matrix$VJ_sequence_nt_trimmed
          if (is.null(is.trimmed) == F) {
            VDJ.matrix$pasted_variants = paste(VDJ.matrix$VDJ_sequence_nt_trimmed, 
                                               VDJ.matrix$VJ_sequence_nt_trimmed, sep = ";")
            message(paste0("Trimmed sequences found for sample ", 
                           sample.names[i]))
            message("Variants are obtained as VDJ_sequence_nt_trimmed;VJ_sequence_nt_trimmed")
          }
          if (is.null(is.trimmed) == T) {
            VDJ.matrix$pasted_variants = paste(VDJ.matrix$VDJ_sequence_nt_raw, 
                                               VDJ.matrix$VJ_sequence_nt_raw, sep = ";")
            message(paste0("Trimmed sequences not found for sample ", 
                           sample.names[i]))
            message("Variants are obtained as VDJ_sequence_nt_raw;VJ_sequence_nt_raw")
          }
          curr_rep_iso <- VDJ.matrix[, c("barcode", 
                                         "clonotype_id", "VDJ_cgene", "pasted_variants")]
          names(curr_rep_iso)[3] <- "isotype"
          curr_rep_iso$isotype <- substr(curr_rep_iso$isotype, 
                                         1, 4)
          if (color.by != "isotype") {
            curr_rep_iso$colors <- as.character(VDJ.matrix[, 
                                                           color.by])
            curr_rep_iso$colors[which(is.na(curr_rep_iso$colors))] <- "None"
            if (inherits(VDJ.matrix[, color.by], "factor")) {
              curr_rep_iso$colors <- ordered(as.factor(curr_rep_iso$colors), 
                                             levels = c(levels(VDJ.matrix[, color.by]), 
                                                        "None"))
            }
            else {
              curr_rep_iso$colors <- ordered(as.factor(curr_rep_iso$colors), 
                                             levels = c(as.character(unique(VDJ.matrix[, 
                                                                                       color.by])), "None"))
            }
          }
          if (treat.incomplete.clones == "exclude") {
            clones_to_del <- c()
            for (k in 1:nrow(clono_freq)) {
              if (all(curr_rep_iso[which(curr_rep_iso$clonotype_id == 
                                         clono_freq[k, 1]), "isotype"] == "")) {
                clones_to_del <- append(clones_to_del, 
                                        k)
              }
            }
            if (length(clones_to_del > 0)) {
              s_id_num = clones_to_del
              s_id_num = sprintf("%05d", s_id_num)
              s_id_num = paste0("clonotype", s_id_num)
              clones_to_del <- s_id_num
              s_id = NULL
              s_id_num = NULL
              s_id = curr_rep_iso$clonotype_id
              s_id_num = readr::parse_number(s_id)
              s_id_num = sprintf("%05d", s_id_num)
              s_id_num = paste0("clonotype", s_id_num)
              curr_rep_iso$clonotype_id <- s_id_num
              s_id = NULL
              s_id_num = NULL
              curr_rep_iso$match = curr_rep_iso$clonotype_id %in% 
                clones_to_del
              curr_rep_iso <- subset(curr_rep_iso, subset = match == 
                                       F)
              curr_rep_iso$match = NULL
              curr_rep_iso = curr_rep_iso[order(curr_rep_iso$clonotype_id), 
              ]
              curr_rep_iso$clonotype_id <- paste0("old", 
                                                  curr_rep_iso$clonotype_id)
              rank_raw <- as.data.frame(curr_rep_iso %>% 
                                          dplyr::group_by(clonotype_id) %>% dplyr::summarise(n = n()))
              rank_raw <- rank_raw[order(rank_raw$n, 
                                         rank_raw$clonotype_id, decreasing = T), 
              ]
              rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
              rank_raw$new_clonotype = paste0("newclonotype", 
                                              sprintf("%05d", rank_raw$new_clonotype))
              for (z in 1:nrow(rank_raw)) {
                curr_rep_iso$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                                  rank_raw$new_clonotype[z], curr_rep_iso$clonotype_id)
              }
              curr_rep_iso <- curr_rep_iso[order(curr_rep_iso$clonotype_id), 
              ]
            }
          }
          s_id = curr_rep_iso$clonotype_id
          s_id_num = readr::parse_number(s_id)
          s_id_num = sprintf("%05d", s_id_num)
          s_id_num = paste0("clonotype", s_id_num)
          curr_rep_iso$clonotype_id <- s_id_num
          s_id = NULL
          s_id_num = NULL
          options(dplyr.summarise.inform = FALSE)
          if (color.by == "isotype") {
            variant_df <- curr_rep_iso %>% dplyr::group_by(pasted_variants, 
                                                           clonotype_id, isotype) %>% dplyr::summarise(n = n())
          }
          if (color.by != "isotype") {
            variant_df <- curr_rep_iso %>% dplyr::group_by(pasted_variants, 
                                                           clonotype_id, colors) %>% dplyr::summarise(n = n())
          }
          variant_df = variant_df[order(variant_df$clonotype_id, 
                                        variant_df$n), ]
          max_row = max(which(variant_df$clonotype_id == 
                                paste0("clonotype", sprintf("%05d", clones))))
          variant_df = variant_df[1:max_row, ]
          s_id = variant_df$clonotype_id
          s_id_num = readr::parse_number(s_id)
          variant_df$clonotype_id <- s_id_num
          s_id = NULL
          s_id_num = NULL
          if (treat.incomplete.cells == "proportional" & 
              stringr::str_detect(paste0(variant_df$isotype, 
                                         collapse = ";"), pattern = "IG") == T) {
            for (clone_number in 1:clones) {
              subset_clone_variant_with_isotype = variant_df[which(variant_df$clonotype_id == 
                                                                     clone_number & stringr::str_detect(variant_df$isotype, 
                                                                                                        pattern = "IG") == T), ]
              subset_clone_variant = variant_df[which(variant_df$clonotype_id == 
                                                        clone_number), ]
              variants <- c()
              if ((nrow(subset_clone_variant_with_isotype) < 
                   nrow(subset_clone_variant)) & (nrow(subset_clone_variant_with_isotype) >= 
                                                  1)) {
                for (z in 1:nrow(subset_clone_variant_with_isotype)) {
                  variants <- c(variants, rep(subset_clone_variant_with_isotype$isotype[z], 
                                              subset_clone_variant_with_isotype$n[z]))
                }
                proportional_variants = c()
                for (variantes in 1:length(which(variant_df$clonotype_id == 
                                                 clone_number & stringr::str_detect(variant_df$isotype, 
                                                                                    pattern = "IG") == F))) {
                  proportional_variants = c(proportional_variants, 
                                            sample(variants, 1))
                }
                variant_df[which(variant_df$clonotype_id == 
                                   clone_number & stringr::str_detect(variant_df$isotype, 
                                                                      pattern = "IG") == F), ]$isotype = proportional_variants
              }
              else if ((nrow(subset_clone_variant_with_isotype) < 
                        nrow(subset_clone_variant)) & (nrow(subset_clone_variant_with_isotype) == 
                                                       0)) {
                proportional_variants = c()
                proportional_variants = c(proportional_variants, 
                                          rep("Unknown", nrow(subset_clone_variant)))
                variant_df[which(variant_df$clonotype_id == 
                                   clone_number & stringr::str_detect(variant_df$isotype, 
                                                                      pattern = "IG") == F), ]$isotype = proportional_variants
              }
            }
          }
          if (treat.incomplete.cells == "exclude" & 
              color.by == "isotype") {
            variant_df$match_ex_crit <- variant_df$isotype == 
              ""
            variant_df <- subset(variant_df, subset = match_ex_crit == 
                                   F)
            variant_df$match_ex_crit = NULL
            variant_df$original_clonotypes = variant_df$clonotype_id
            variant_df$clonotype_id = sprintf("%05d", 
                                              variant_df$clonotype_id)
            variant_df$clonotype_id = paste0("old", 
                                             variant_df$clonotype_id)
            rank_raw <- as.data.frame(variant_df %>% 
                                        dplyr::group_by(clonotype_id) %>% dplyr::summarise(n = sum(n)))
            rank_raw <- rank_raw[order(rank_raw$n, rank_raw$clonotype_id, 
                                       decreasing = T), ]
            rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
            rank_raw$new_clonotype = sprintf("%05d", 
                                             rank_raw$new_clonotype)
            rank_raw$new_clonotype = paste0("new", rank_raw$new_clonotype)
            for (z in 1:nrow(rank_raw)) {
              variant_df$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                              rank_raw$new_clonotype[z], variant_df$clonotype_id)
            }
            variant_df$clonotype_id <- readr::parse_number(variant_df$clonotype_id)
            variant_df <- variant_df[order(variant_df$clonotype_id), 
            ]
          }
          if (color.by == "isotype" & isotypes.to.plot[1] != 
              "all") {
            if (!any(isotypes.to.plot %in% unique(variant_df$isotype))) {
              stop("isotype.to.plot input not found in dataframe. Please check if the isotype is spelled correctly")
            }
            variant_df$match_ex_crit <- variant_df$isotype %in% 
              isotypes.to.plot
            variant_df <- subset(variant_df, subset = match_ex_crit == 
                                   T)
            variant_df$match_ex_crit = NULL
            variant_df$original_clonotypes = variant_df$clonotype_id
            variant_df$clonotype_id = sprintf("%05d", 
                                              variant_df$clonotype_id)
            variant_df$clonotype_id = paste0("old", 
                                             variant_df$clonotype_id)
            rank_raw <- as.data.frame(variant_df %>% 
                                        dplyr::group_by(clonotype_id) %>% dplyr::summarise(n = sum(n)))
            rank_raw <- rank_raw[order(rank_raw$n, rank_raw$clonotype_id, 
                                       decreasing = T), ]
            rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
            rank_raw$new_clonotype = sprintf("%05d", 
                                             rank_raw$new_clonotype)
            rank_raw$new_clonotype = paste0("new", rank_raw$new_clonotype)
            for (z in 1:nrow(rank_raw)) {
              variant_df$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                              rank_raw$new_clonotype[z], variant_df$clonotype_id)
            }
            variant_df$clonotype_id <- readr::parse_number(variant_df$clonotype_id)
            variant_df <- variant_df[order(variant_df$clonotype_id), 
            ]
            message("New ranking based only on selected isotypes: ")
            message(paste0(unique(variant_df$clonotype_id), 
                           ", "))
          }
          variant_df$variant = 1:length(variant_df$clonotype_id)
          variant_df$variant = as.character(variant_df$variant)
          variant_df_list[[i]] = variant_df
          if (color.by == "isotype") {
            output_plot[[i]] <- ggplot2::ggplot(variant_df_list[[i]], 
                                                ggplot2::aes(x = clonotype_id, y = n, 
                                                             fill = isotype, color = variant)) + 
              ggplot2::geom_bar(stat = "identity", width = 0.6, 
                                color = "white") + ggplot2::theme_bw() + 
              ggplot2::scale_fill_manual("Isotype", 
                                         values = c(IGHG = "green4", IGHM = "black", 
                                                    IGHA = "red3", IGHD = "blue", IGHE = "purple", 
                                                    Unknown = "gray")) + ggplot2::theme_classic() + 
              ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                                  0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                        x = "Clonal rank", y = "Number of cells")
          }
          else {
            output_plot[[i]] <- ggplot2::ggplot(variant_df_list[[i]], 
                                                ggplot2::aes(x = clonotype_id, y = n, 
                                                             fill = colors, color = variant)) + ggplot2::geom_bar(stat = "identity", 
                                                                                                                  width = 0.6, color = "white") + ggplot2::scale_fill_manual(values = grDevices::rainbow(n = length(unique(clones_per_isotype_all[[i]]$colors)))) + 
              ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + 
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                                  0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                        x = "Clonal rank", y = "Number of cells")
          }
        }
        clones_per_isotype_all = variant_df_list
        names(clones_per_isotype_all) <- sample.names
        names(output_plot) <- sample.names
        return(list(output_plot, clones_per_isotype_all))
      }
    }
    if (subtypes == TRUE & color.by == "isotype") {
      if (variant.plot == FALSE) {
        for (i in 1:length(sample.names)) {
          VDJ.matrix <- subset(VDJ.matrix.all, VDJ.matrix.all[, 
                                                              group.by] == sample.names[i])
          clono_freq <- as.data.frame(table(VDJ.matrix$clonotype_id))
          clono_freq <- clono_freq[order(clono_freq$Freq, 
                                         decreasing = T), ]
          curr_rep_iso <- VDJ.matrix[, c("barcode", 
                                         "clonotype_id", "VDJ_cgene", "VDJ_cdr3s_aa", 
                                         "VJ_cdr3s_aa")]
          curr_rep_iso$VDJ_cgene[is.null(curr_rep_iso$VDJ_cgene)] <- "none"
          names(curr_rep_iso)[3] <- "isotype"
          curr_rep_iso$isotype <- stringr::str_split(curr_rep_iso$isotype, 
                                                     ";", simplify = T)[, 1]
          if (treat.incomplete.clones == "exclude") {
            clones_to_del <- c()
            for (k in 1:nrow(clono_freq)) {
              if (all(curr_rep_iso[which(curr_rep_iso$clonotype_id == 
                                         clono_freq[k, 1]), "isotype"] == "")) {
                clones_to_del <- append(clones_to_del, 
                                        k)
              }
            }
            if (length(clones_to_del > 0)) {
              clono_freq <- clono_freq[-clones_to_del, 
              ]
            }
            clono_freq <- clono_freq[order(clono_freq$Freq, 
                                           decreasing = T), ]
          }
          clones_per_isotype <- list()
          j <- 1
          for (j in 1:clones) {
            curr_clone <- curr_rep_iso[which(curr_rep_iso$clonotype_id == 
                                               clono_freq[j, 1]), ]
            if (nrow(curr_clone) > 0) {
              curr_clone$isotype <- stringr::str_split(curr_clone$isotype, 
                                                       ";", simplify = T)[, 1]
              if (treat.incomplete.cells == "proportional" & 
                  stringr::str_detect(paste0(curr_clone$isotype, 
                                             collapse = ";"), pattern = "IG") == 
                  T) {
                props <- table(curr_clone$isotype[which(stringr::str_detect(curr_clone$isotype, 
                                                                            pattern = "IG") == T)])
                n_total <- nrow(curr_clone)
                props <- round(props/sum(props) * n_total, 
                               0)
                if (sum(props) > n_total) {
                  props[which.max(props)] <- props[which.max(props)] - 
                    1
                }
                else if (sum(props) < n_total) {
                  props[which.min(props)] <- props[which.min(props)] + 
                    1
                }
                curr_clone$isotype <- rep.int(names(props), 
                                              props)
              }
              else if (treat.incomplete.cells == "proportional" & 
                       stringr::str_detect(paste0(curr_clone$isotype, 
                                                  collapse = ";"), pattern = "IG") == 
                       F) {
                curr_clone$isotype <- "None"
              }
              clones_per_isotype[[j]] <- data.frame(Counts = rep(0, 
                                                                 14), Color = rep("", 14), Isotype = rep("", 
                                                                                                         14), ClonalRank = rep("", 14), clonotype_id = rep(curr_clone$clonotype_id[1], 
                                                                                                                                                           14), VDJ_cdr3s_aa = rep(curr_clone$VDJ_cdr3s_aa[which(curr_clone$VDJ_cdr3s_aa != 
                                                                                                                                                                                                                   "")][1], 14), VJ_cdr3s_aa = rep(curr_clone$VJ_cdr3s_aa[which(curr_clone$VJ_cdr3s_aa != 
                                                                                                                                                                                                                                                                                  "")][1], 14), barcode = rep(paste0(curr_clone$barcode, 
                                                                                                                                                                                                                                                                                                                     collapse = ";"), 14))
              clones_per_isotype[[j]]$Counts[1] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHG1"))
              clones_per_isotype[[j]]$Counts[2] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHG2"))
              clones_per_isotype[[j]]$Counts[3] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHG2A"))
              clones_per_isotype[[j]]$Counts[4] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHG2B"))
              clones_per_isotype[[j]]$Counts[5] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHG2C"))
              clones_per_isotype[[j]]$Counts[6] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHG3"))
              clones_per_isotype[[j]]$Counts[7] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHG4"))
              clones_per_isotype[[j]]$Counts[8] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHM"))
              clones_per_isotype[[j]]$Counts[9] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                          "IGHA"))
              clones_per_isotype[[j]]$Counts[10] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                           "IGHA1"))
              clones_per_isotype[[j]]$Counts[11] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                           "IGHA2"))
              clones_per_isotype[[j]]$Counts[12] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                           "IGHD"))
              clones_per_isotype[[j]]$Counts[13] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                           "IGHE"))
              clones_per_isotype[[j]]$Counts[14] <- sum(stringr::str_count(curr_clone$isotype, 
                                                                           "None"))
              clones_per_isotype[[j]]$ClonalRank <- j
              if (color.by == "isotype") {
                clones_per_isotype[[j]]$Isotype <- c("IGHG1", 
                                                     "IGHG2", "IGHG2a", "IGHG2b", "IGHG2c", 
                                                     "IGHG3", "IGHG4", "IGHM", "IGHA", 
                                                     "IGHA1", "IGHA2", "IGHD", "IGHE", 
                                                     "Unknown")
                if (species == "Human") {
                  clones_per_isotype[[j]] <- clones_per_isotype[[j]][-c(3, 
                                                                        4, 5, 9), ]
                }
                if (species == "Mouse") {
                  clones_per_isotype[[j]] <- clones_per_isotype[[j]][-c(2, 
                                                                        7, 10, 11), ]
                }
              }
              else {
                clones_per_isotype[[j]]$Isotype <- names(which.max(table(curr_clone$colors)))
              }
            }
          }
          clones_per_isotype_all[[i]] <- do.call("rbind", 
                                                 clones_per_isotype)
          if (treat.incomplete.cells == "exclude") {
            rank_raw <- as.data.frame(clones_per_isotype_all[[i]] %>% 
                                        dplyr::group_by(ClonalRank) %>% dplyr::summarise(sum_counts = sum(Counts)) %>% 
                                        dplyr::arrange(dplyr::desc(sum_counts)) %>% 
                                        dplyr::mutate(rank = 1:length(unique(ClonalRank))))
            clones_per_isotype_all[[i]]$ClonalRank_2 <- 0
            for (l in 1:nrow(rank_raw)) {
              clones_per_isotype_all[[i]]$ClonalRank_2[which(clones_per_isotype_all[[i]]$ClonalRank == 
                                                               rank_raw$ClonalRank[l])] <- rank_raw$rank[l]
            }
            clones_per_isotype_all[[i]]$ClonalRank <- clones_per_isotype_all[[i]]$ClonalRank_2
            message("New ranking based only on present HC chains: ")
            message(unique(clones_per_isotype_all[[i]]$clonotype_id))
          }
          if (color.by == "isotype" & isotypes.to.plot[1] != 
              "all") {
            if (!any(isotypes.to.plot %in% unique(clones_per_isotype_all[[i]]$Isotype))) {
              stop("isotype.to.plot input not found in dataframe. Please check if the isotype is spelled correctly")
            }
            to_del <- c()
            for (k in 1:length(unique(clones_per_isotype_all[[i]]$ClonalRank))) {
              if (sum(clones_per_isotype_all[[i]]$Counts[which(clones_per_isotype_all[[i]]$ClonalRank == 
                                                               unique(clones_per_isotype_all[[i]]$ClonalRank)[k] & 
                                                               clones_per_isotype_all[[i]]$Isotype %in% 
                                                               isotypes.to.plot)]) == 0) {
                to_del <- append(to_del, unique(clones_per_isotype_all[[i]]$ClonalRank)[k])
              }
            }
            clones_per_isotype_all[[i]] <- subset(clones_per_isotype_all[[i]], 
                                                  !ClonalRank %in% to_del)
            rank_raw <- as.data.frame(clones_per_isotype_all[[i]] %>% 
                                        dplyr::group_by(ClonalRank) %>% dplyr::summarise(sum_counts = sum(Counts)) %>% 
                                        plyr::arrange(dplyr::desc(sum_counts)) %>% 
                                        dplyr::mutate(rank = 1:length(unique(ClonalRank))))
            clones_per_isotype_all[[i]]$ClonalRank_2 <- 0
            for (l in 1:nrow(rank_raw)) {
              clones_per_isotype_all[[i]]$ClonalRank_2[which(clones_per_isotype_all[[i]]$ClonalRank == 
                                                               rank_raw$ClonalRank[l])] <- rank_raw$rank[l]
            }
            clones_per_isotype_all[[i]]$ClonalRank <- clones_per_isotype_all[[i]]$ClonalRank_2
            message("New ranking based only on selected isotypes: ")
            message(unique(clones_per_isotype_all[[i]]$clonotype_id))
          }
          if (species == "Human") {
            
            
            output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], 
                                                ggplot2::aes(fill = Isotype, y = Counts, 
                                                             x = ClonalRank)) + ggplot2::geom_bar(stat = "identity", 
                                                                                                  width = 0.6, color = "black") + ggplot2::theme_bw() + 
              ggplot2::scale_fill_manual("Isotype", 
                                         values = c(IGHG1 = "green", IGHG2 = "green3", 
                                                    IGHG3 = "green4", IGHG4 = "darkgreen", 
                                                    IGHM = "black", IGHA1 = "red", IGHA2 = "red4", 
                                                    IGHD = "blue", IGHE = "purple", Unknown = "gray")) + 
              ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + 
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::ylab("Number of cells") + 
              ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                    0.5))
          }
          if (species == "Mouse") {
            output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], 
                                                ggplot2::aes(fill = Isotype, y = Counts, 
                                                             x = ClonalRank)) + ggplot2::geom_bar(stat = "identity", 
                                                                                                  width = 0.6, color = "black") + ggplot2::theme_bw() + 
              ggplot2::scale_fill_manual("Isotype", 
                                         values = c(IGHG1 = "lightgreen", IGHG2a = "green", 
                                                    IGHG2b = "green3", IGHG2c = "green4", 
                                                    IGHG3 = "darkgreen", IGHM = "black", 
                                                    IGHA = "red", IGHD = "blue", IGHE = "purple", 
                                                    Unknown = "gray")) + ggplot2::theme_classic() + 
              ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::ylab("Number of cells") + 
              ggplot2::xlab("Clonal rank") + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                    0.5))
          }
        }
        names(clones_per_isotype_all) <- sample.names
        return(list(output_plot, clones_per_isotype_all))
      }
      if (variant.plot == TRUE) {
        for (i in 1:length(sample.names)) {
          message(paste0("Starting sample ", i, "/", 
                         length(sample.names)))
          VDJ.matrix <- subset(VDJ.matrix.all, VDJ.matrix.all[, 
                                                              group.by] == sample.names[i])
          s_id = VDJ.matrix$clonotype_id
          s_id_num = readr::parse_number(s_id)
          s_id_num = sprintf("%05d", s_id_num)
          s_id_num = paste0("clonotype", s_id_num)
          VDJ.matrix$clonotype_id <- s_id_num
          s_id = NULL
          s_id_num = NULL
          VDJ.matrix$clonotype_id = paste0("old", VDJ.matrix$clonotype_id)
          rank_raw <- as.data.frame(VDJ.matrix %>% dplyr::group_by(clonotype_id) %>% 
                                      dplyr::summarise(n = n()))
          rank_raw <- rank_raw[order(rank_raw$n, rank_raw$clonotype_id, 
                                     decreasing = T), ]
          rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
          rank_raw$new_clonotype = sprintf("%05d", rank_raw$new_clonotype)
          rank_raw$new_clonotype = paste0("new", rank_raw$new_clonotype)
          for (z in 1:nrow(rank_raw)) {
            VDJ.matrix$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                            rank_raw$new_clonotype[z], VDJ.matrix$clonotype_id)
          }
          VDJ.matrix$clonotype_id <- paste0("clonotype", 
                                            sprintf("%05d", readr::parse_number(VDJ.matrix$clonotype_id)))
          clono_freq <- as.data.frame(table(VDJ.matrix$clonotype_id))
          clono_freq <- clono_freq[order(clono_freq$Freq, 
                                         decreasing = T), ]
          is.trimmed = VDJ.matrix$VJ_sequence_nt_trimmed
          if (is.null(is.trimmed) == F) {
            VDJ.matrix$pasted_variants = paste(VDJ.matrix$VDJ_sequence_nt_trimmed, 
                                               VDJ.matrix$VJ_sequence_nt_trimmed, sep = ";")
            message(paste0("Trimmed sequences found for sample ", 
                           sample.names[i]))
            message("Variants are obtained as VDJ_sequence_nt_trimmed;VJ_sequence_nt_trimmed")
          }
          if (is.null(is.trimmed) == T) {
            VDJ.matrix$pasted_variants = paste(VDJ.matrix$VDJ_sequence_nt_raw, 
                                               VDJ.matrix$VJ_sequence_nt_raw, sep = ";")
            message(paste0("Trimmed sequences not found for sample ", 
                           sample.names[i]))
            message("Variants are obtained as VDJ_sequence_nt_raw;VJ_sequence_nt_raw")
          }
          curr_rep_iso <- VDJ.matrix[, c("barcode", 
                                         "clonotype_id", "VDJ_cgene", "pasted_variants")]
          names(curr_rep_iso)[3] <- "isotype"
          curr_rep_iso$isotype <- stringr::str_split(curr_rep_iso$isotype, 
                                                     ";", simplify = T)[, 1]
          if (treat.incomplete.clones == "exclude") {
            clones_to_del <- c()
            for (k in 1:nrow(clono_freq)) {
              if (all(curr_rep_iso[which(curr_rep_iso$clonotype_id == 
                                         clono_freq[k, 1]), "isotype"] == "")) {
                clones_to_del <- append(clones_to_del, 
                                        k)
              }
            }
            if (length(clones_to_del > 0)) {
              s_id_num = clones_to_del
              s_id_num = sprintf("%05d", s_id_num)
              s_id_num = paste0("clonotype", s_id_num)
              clones_to_del <- s_id_num
              s_id = NULL
              s_id_num = NULL
              s_id = curr_rep_iso$clonotype_id
              s_id_num = readr::parse_number(s_id)
              s_id_num = sprintf("%05d", s_id_num)
              s_id_num = paste0("clonotype", s_id_num)
              curr_rep_iso$clonotype_id <- s_id_num
              s_id = NULL
              s_id_num = NULL
              curr_rep_iso$match = curr_rep_iso$clonotype_id %in% 
                clones_to_del
              curr_rep_iso <- subset(curr_rep_iso, subset = match == 
                                       F)
              curr_rep_iso$match = NULL
              curr_rep_iso = curr_rep_iso[order(curr_rep_iso$clonotype_id), 
              ]
              curr_rep_iso$clonotype_id <- paste0("old", 
                                                  curr_rep_iso$clonotype_id)
              rank_raw <- as.data.frame(curr_rep_iso %>% 
                                          dplyr::group_by(clonotype_id) %>% dplyr::summarise(n = n()))
              rank_raw <- rank_raw[order(rank_raw$n, 
                                         rank_raw$clonotype_id, decreasing = T), 
              ]
              rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
              rank_raw$new_clonotype = paste0("newclonotype", 
                                              sprintf("%05d", rank_raw$new_clonotype))
              for (z in 1:nrow(rank_raw)) {
                curr_rep_iso$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                                  rank_raw$new_clonotype[z], curr_rep_iso$clonotype_id)
              }
              curr_rep_iso <- curr_rep_iso[order(curr_rep_iso$clonotype_id), 
              ]
            }
          }
          s_id = curr_rep_iso$clonotype_id
          s_id_num = readr::parse_number(s_id)
          s_id_num = sprintf("%05d", s_id_num)
          s_id_num = paste0("clonotype", s_id_num)
          curr_rep_iso$clonotype_id <- s_id_num
          s_id = NULL
          s_id_num = NULL
          options(dplyr.summarise.inform = FALSE)
          variant_df <- curr_rep_iso %>% dplyr::group_by(pasted_variants, 
                                                         clonotype_id, isotype) %>% dplyr::summarise(n = n())
          variant_df = variant_df[order(variant_df$clonotype_id, 
                                        variant_df$n), ]
          max_row = max(which(variant_df$clonotype_id == 
                                paste0("clonotype", sprintf("%05d", clones))))
          variant_df = variant_df[1:max_row, ]
          s_id = variant_df$clonotype_id
          s_id_num = readr::parse_number(s_id)
          variant_df$clonotype_id <- s_id_num
          s_id = NULL
          s_id_num = NULL
          if (treat.incomplete.cells == "proportional" & 
              stringr::str_detect(paste0(variant_df$isotype, 
                                         collapse = ";"), pattern = "IG") == T) {
            for (clone_number in 1:clones) {
              subset_clone_variant_with_isotype = variant_df[which(variant_df$clonotype_id == 
                                                                     clone_number & stringr::str_detect(variant_df$isotype, 
                                                                                                        pattern = "IG") == T), ]
              subset_clone_variant = variant_df[which(variant_df$clonotype_id == 
                                                        clone_number), ]
              variants <- c()
              if ((nrow(subset_clone_variant_with_isotype) < 
                   nrow(subset_clone_variant)) & (nrow(subset_clone_variant_with_isotype) >= 
                                                  1)) {
                for (z in 1:nrow(subset_clone_variant_with_isotype)) {
                  variants <- c(variants, rep(subset_clone_variant_with_isotype$isotype[z], 
                                              subset_clone_variant_with_isotype$n[z]))
                }
                proportional_variants = c()
                for (variantes in 1:length(which(variant_df$clonotype_id == 
                                                 clone_number & stringr::str_detect(variant_df$isotype, 
                                                                                    pattern = "IG") == F))) {
                  proportional_variants = c(proportional_variants, 
                                            sample(variants, 1))
                }
                variant_df[which(variant_df$clonotype_id == 
                                   clone_number & stringr::str_detect(variant_df$isotype, 
                                                                      pattern = "IG") == F), ]$isotype = proportional_variants
              }
              else if ((nrow(subset_clone_variant_with_isotype) < 
                        nrow(subset_clone_variant)) & (nrow(subset_clone_variant_with_isotype) == 
                                                       0)) {
                proportional_variants = c()
                proportional_variants = c(proportional_variants, 
                                          rep("Unknown", nrow(subset_clone_variant)))
                variant_df[which(variant_df$clonotype_id == 
                                   clone_number & stringr::str_detect(variant_df$isotype, 
                                                                      pattern = "IG") == F), ]$isotype = proportional_variants
              }
            }
          }
          if (treat.incomplete.cells == "exclude" & 
              color.by == "isotype") {
            variant_df$match_ex_crit <- variant_df$isotype == 
              ""
            variant_df <- subset(variant_df, subset = match_ex_crit == 
                                   F)
            variant_df$match_ex_crit = NULL
            variant_df$original_clonotypes = variant_df$clonotype_id
            variant_df$clonotype_id = sprintf("%05d", 
                                              variant_df$clonotype_id)
            variant_df$clonotype_id = paste0("old", 
                                             variant_df$clonotype_id)
            rank_raw <- as.data.frame(variant_df %>% 
                                        dplyr::group_by(clonotype_id) %>% dplyr::summarise(n = sum(n)))
            rank_raw <- rank_raw[order(rank_raw$n, rank_raw$clonotype_id, 
                                       decreasing = T), ]
            rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
            rank_raw$new_clonotype = sprintf("%05d", 
                                             rank_raw$new_clonotype)
            rank_raw$new_clonotype = paste0("new", rank_raw$new_clonotype)
            for (z in 1:nrow(rank_raw)) {
              variant_df$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                              rank_raw$new_clonotype[z], variant_df$clonotype_id)
            }
            variant_df$clonotype_id <- readr::parse_number(variant_df$clonotype_id)
            variant_df <- variant_df[order(variant_df$clonotype_id), 
            ]
          }
          if (color.by == "isotype" & isotypes.to.plot[1] != 
              "all") {
            if (!any(isotypes.to.plot %in% unique(variant_df$isotype))) {
              stop("isotype.to.plot input not found in dataframe. Please check if the isotype is spelled correctly")
            }
            variant_df$match_ex_crit <- variant_df$isotype %in% 
              isotypes.to.plot
            variant_df <- subset(variant_df, subset = match_ex_crit == 
                                   T)
            variant_df$match_ex_crit = NULL
            variant_df$original_clonotypes = variant_df$clonotype_id
            variant_df$clonotype_id = sprintf("%05d", 
                                              variant_df$clonotype_id)
            variant_df$clonotype_id = paste0("old", 
                                             variant_df$clonotype_id)
            rank_raw <- as.data.frame(variant_df %>% 
                                        dplyr::group_by(clonotype_id) %>% dplyr::summarise(n = sum(n)))
            rank_raw <- rank_raw[order(rank_raw$n, rank_raw$clonotype_id, 
                                       decreasing = T), ]
            rank_raw$new_clonotype <- 1:length(rank_raw$clonotype_id)
            rank_raw$new_clonotype = sprintf("%05d", 
                                             rank_raw$new_clonotype)
            rank_raw$new_clonotype = paste0("new", rank_raw$new_clonotype)
            for (z in 1:nrow(rank_raw)) {
              variant_df$clonotype_id <- gsub(rank_raw$clonotype_id[z], 
                                              rank_raw$new_clonotype[z], variant_df$clonotype_id)
            }
            variant_df$clonotype_id <- readr::parse_number(variant_df$clonotype_id)
            variant_df <- variant_df[order(variant_df$clonotype_id), 
            ]
            message("New ranking based only on selected isotypes: ")
            message(paste0(unique(variant_df$clonotype_id), 
                           ", "))
          }
          variant_df$variant = 1:length(variant_df$clonotype_id)
          variant_df$variant = as.character(variant_df$variant)
          variant_df_list[[i]] = variant_df
          if (species == "Human") {
            output_plot[[i]] <- ggplot2::ggplot(variant_df_list[[i]], 
                                                ggplot2::aes(x = clonotype_id, y = n, 
                                                             fill = isotype, color = variant)) + 
              ggplot2::geom_bar(stat = "identity", width = 0.6, 
                                color = "white") + ggplot2::theme_bw() + 
              ggplot2::scale_fill_manual("Isotype", 
                                         values = c(IGHG1 = "green", IGHG2 = "green3", 
                                                    IGHG3 = "green4", IGHG4 = "darkgreen", 
                                                    IGHM = "black", IGHA1 = "red", IGHA2 = "red4", 
                                                    IGHD = "blue", IGHE = "purple", Unknown = "gray")) + 
              ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + 
              ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                                  0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                        x = "Clonal rank", y = "Number of cells")
          }
          if (species == "Mouse") {
            output_plot[[i]] <- ggplot2::ggplot(variant_df_list[[i]], 
                                                ggplot2::aes(x = clonotype_id, y = n, 
                                                             fill = isotype, color = variant)) + 
              ggplot2::geom_bar(stat = "identity", width = 0.6, 
                                color = "white") + ggplot2::theme_bw() + 
              ggplot2::scale_fill_manual("Isotype", 
                                         values = c(IGHG1 = "lightgreen", IGHG2A = "green", 
                                                    IGHG2B = "green3", IGHG2C = "green4", 
                                                    IGHG3 = "darkgreen", IGHM = "black", 
                                                    IGHA = "red", IGHD = "blue", IGHE = "purple", 
                                                    Unknown = "gray")) + ggplot2::theme_classic() + 
              ggplot2::ggtitle(paste0(i)) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
              ggplot2::scale_y_continuous(expand = c(0, 
                                                     0)) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                                  0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                        x = "Clonal rank", y = "Number of cells")
          }
        }
        clones_per_isotype_all = variant_df_list
        names(clones_per_isotype_all) <- sample.names
        names(output_plot) <- sample.names
        return(list(output_plot, clones_per_isotype_all))
      }
    }
  }  else if (celltype == "Tcells") {
    if (!inherits(VDJ.matrix, "data.frame")) {
      stop("Please provide a VDJ matrix dataframe (VDJ_GEX_matrix_output[[1]])")
    }
    VDJ.matrix[, group.by] <- as.character(VDJ.matrix[, 
                                                      group.by])
    sample.names <- unique(VDJ.matrix[, group.by])
    sample.names[is.na(sample.names)] <- "NONE"
    VDJ.matrix[is.na(VDJ.matrix[, group.by]), group.by] <- "NONE"
    VDJ.matrix.all <- VDJ.matrix
    if (clones < 1 | clones > 50) {
      stop("Number of clones must be an integer value between 1 and 50")
    }
    VDJ_per_clone_output_all <- list()
    clones_per_isotype <- list()
    clones_per_isotype_all <- list()
    output_plot <- list()
    for (i in 1:length(sample.names)) {
      VDJ.matrix <- subset(VDJ.matrix.all, VDJ.matrix.all[, 
                                                          group.by] == sample.names[i])
      clono_freq <- as.data.frame(table(VDJ.matrix$clonotype_id))
      clono_freq <- clono_freq[order(clono_freq$Freq, 
                                     decreasing = T), ]
      curr_rep_iso <- VDJ.matrix[, c("barcode", "clonotype_id", 
                                     "VDJ_cdr3s_aa", "VJ_cdr3s_aa")]
      if (color.by[1] != "isotype") {
        curr_rep_iso$colors <- as.character(VDJ.matrix[, 
                                                       color.by])
        curr_rep_iso$colors[which(is.na(curr_rep_iso$colors))] <- "None"
        if (inherits(VDJ.matrix[, color.by], "factor")) {
          curr_rep_iso$colors <- ordered(as.factor(curr_rep_iso$colors), 
                                         levels = c(levels(VDJ.matrix[, color.by]), 
                                                    "None"))
        }
        else {
          curr_rep_iso$colors <- ordered(as.factor(curr_rep_iso$colors), 
                                         levels = c(as.character(unique(VDJ.matrix[, 
                                                                                   color.by])), "None"))
        }
      }
      else {
        curr_rep_iso$colors <- "none"
      }
      curr_rep_iso$isotype <- "none"
      clones_per_isotype <- list()
      for (j in 1:clones) {
        curr_clone <- curr_rep_iso[which(curr_rep_iso$clonotype_id == 
                                           clono_freq[j, 1]), ]
        if (nrow(curr_clone) > 0) {
          color_cur_clone <- unique(curr_clone$colors)
          n_color_cur_clone <- length(unique(curr_clone$colors))
          clones_per_isotype[[j]] <- data.frame(Counts = rep(0, 
                                                             n_color_cur_clone), Color = color_cur_clone, 
                                                ClonalRank = rep("", n_color_cur_clone), 
                                                clonotype_id = rep(curr_clone$clonotype_id[1], 
                                                                   n_color_cur_clone), VDJ_cdr3s_aa = rep(curr_clone$VDJ_cdr3s_aa[which(curr_clone$VDJ_cdr3s_aa != 
                                                                                                                                          "")][1], n_color_cur_clone), VJ_cdr3s_aa = rep(curr_clone$VJ_cdr3s_aa[which(curr_clone$VJ_cdr3s_aa != 
                                                                                                                                                                                                                        "")][1], n_color_cur_clone), barcode = rep(paste0(curr_clone$barcode, 
                                                                                                                                                                                                                                                                          collapse = ";"), n_color_cur_clone))
          if (color.by[1] == "isotype") {
            clones_per_isotype[[j]]$Color <- c("None")
            clones_per_isotype[[j]]$Counts <- nrow(curr_clone)
          }
          else {
            for (k in 1:nrow(clones_per_isotype[[j]])) {
              clones_per_isotype[[j]]$Counts[k] <- stringr::str_count(paste0("/", 
                                                                             paste0(curr_clone$colors, collapse = "/ /"), 
                                                                             "/"), pattern = paste0("/", as.character(clones_per_isotype[[j]]$Color[k]), 
                                                                                                    "/"))
            }
          }
          clones_per_isotype[[j]]$ClonalRank <- j
        }
      }
      clones_per_isotype_all[[i]] <- do.call("rbind", 
                                             clones_per_isotype)
      if (color.by[1] == "isotype") {
        output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], 
                                            ggplot2::aes(fill = Color, y = Counts, x = ClonalRank)) + 
          ggplot2::geom_bar(stat = "identity", width = 0.6, 
                            color = "black") + ggplot2::theme_bw() + 
          ggplot2::scale_fill_manual(values = c("gray80")) + 
          ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + 
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
          ggplot2::scale_y_continuous(expand = c(0, 
                                                 0)) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                              0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                    x = "Clonal rank", y = "Number of cells")
      }
      else {
        output_plot[[i]] <- ggplot2::ggplot(clones_per_isotype_all[[i]], 
                                            ggplot2::aes(fill = Color, y = Counts, x = ClonalRank)) + 
          ggplot2::geom_bar(stat = "identity", width = 0.6, 
                            color = "black") + ggplot2::theme_bw() + 
          ggplot2::theme_classic() + ggplot2::ggtitle(paste0(i)) + 
          ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) + 
          ggplot2::scale_y_continuous(expand = c(0, 
                                                 0)) + ggplot2::scale_x_continuous(expand = c(0, 
                                                                                              0.5)) + ggplot2::labs(title = sample.names[[i]], 
                                                                                                                    x = "Clonal rank", y = "Number of cells", 
                                                                                                                    fill = color.by) + ggplot2::scale_fill_manual(values = grDevices::rainbow(n = length(unique(clones_per_isotype_all[[i]]$Color))))
      }
    }
    names(clones_per_isotype_all) <- sample.names
    
    return(list(output_plot, clones_per_isotype_all))
  }
}