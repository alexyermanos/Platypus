#' Simulate repertoire and transcriptome matrix, with igraph tree plot for each clone showing the evolution process. the node in the tree plot are colored with transcriptome state and isotype.
#'
#' @title Simulate immune repertoire and transcriptome data
#'
#' @param initial.size.of.repertoire The initial number of existing cells when the evolution starts. Default is 10.
#' @param species The species of the simulated repertoire, can be "mus" for mouse or "hum" for human. Default is "mus".
#' @param cell.type The cell type for the simulation. "B" or "T"
#' @param cd4.proportion A number between 0 and 1 specifying the proportion of Cd4+ T cells, when cell.type is "T" and transcirptime states data is default. Default is 1, all the cells are Cd4. When user specify transciptome data for T cells, mixture of CD4+ and CD8+ T cells are not applicable.
#' @param duration.of.evolution The maxim time steps for simulation.
#' @param complete.duration TRUE or FALSE. Default is TURE. If TURE, after cell number or clone number reaches the upper limit, the evolution(class switch, mutation, transcriptional state switch) will continue until the duration.of.evolution is complete. If FLASE, the evolution will stop when either cell number or clone number reaches the limit.
#' @param vdj.productive "random": the sequence will be generated from random VDJ recombination, there might be a proportion of unproductive sequences.These VDJ genes were taken from IMGT. When more than one allele was present for a given gene, the first was used. "naive": the VDJ sequence be sampled from a pool of productive sequences obtained by filtering randomly simulated sequences with MIXCR. "vae": the VDJ sequence be sampled from a pool of productive sequences obtained by filtering sequences generated from vae models with MIXCR.
#' @param vdj.model Specifies the model used to simulate V-D-J recombination. Can be either "naive" or "data". "naive" is chain independent and does not differentiate between different species. To rely on the default "experimental" options, this should be "data" and the parameter vdj.insertion.mean should be "default". This will allow for different mean additions for either the VD and JD junctions and will differ depending on species.
#' @param vdj.insertion.mean Integer value describing the mean number of nucleotides to be inserted during simulated V-D-J recombination events. If "default" is entered, the mean will be normally distributed.
#' @param vdj.insertion.stdv Integer value describing the standard deviation corresponding to insertions of V-D-J recombination. No "default" parameter currently supported but will be updated with future experimental data. This should be a number if using a custom distribution for V-D-J recombination events, but can be "default" if using the "naive" vdj.model or the "data", with vdj.insertion.mean set to "default".
#' @param vdj.branch.prob Probability of new VDJ recombination event in each time step.when new VDJ recombination happen, a new cell with a new sequence will be generated. Default is 0.2.
#' @param clonal.selection TRUE or FALSE. If TURE, cells in clones with higher frequency have their division probability proportional to the clonal frequency. If FALSE, clones with higher frequency will have lower probability to expand.
#' @param cell.division.prob Probability of cells to be duplicated in each time step. Default is 0.1. If uneven probability for different clones is needed, the input should be a vector of 2 numeric items, with the first item being the lower bound, the second item being the upper bound of the division rate. The most abundant clone will get the highest division rate, and division rate of other clones will follow arithmetic progression and keep decreasing until the last abundant clone with the lower limit of division rate. If input 3 values, the third value will be the division rate for cells with selected sequences. If a fourth number is given, the division probability of selected sequence will be sampled between the third number and the fourth number.
#' @param sequence.selection.prob Probability of each unique sequence to be selected as expanding sequence.Expanding sequences can have their division rate specified in the third element of cell.division.prob.
#' @param special.v.gene If TRUE, simulation will apply special sequence.selection.prob for heavy and light chain v gene combination specified in dataframe "special_v".
#' @param class.switch.prob Probability matrix of class switching for b cells. The row names of the matrix are the isotypes the cell is switching from, the column names are the isotypes the cell is switching to. All B cells start from IGHM, and switch to one of the other isotypes or remain the same. Default values are in the attaching matrix "class_switch_prob_hum" and "class_switch_prob_mus". The order of isotype in rows and columns should be the same.
#' @param class.switch.selection.dependent If TRUE,class switching will happen when the cell is selected,if the cell has IgM or IgD isotype.
#' @param class.switch.independent If TRUE, class switching will happen randomly at each time step for all cells. If FALSE, random class switching will be switched off.
#' @param SHM.method The mode of SHM speciation events. Options are either: "poisson","data","motif","wrc", and "all". Specifying either "poisson" or "naive" will result in mutations that can occur anywhere in the heavy chain region, with each nucleotide having an equal probability for a mutation event. Specifying "data" focuses mutation events during SHM in the CDR regions (based on IMGT), and there will be an increased probability for transitions (and decreased probability for transversions). Specifying "motif" will cause neighbor dependent mutations based on a mutational matrix from high throughput sequencing data sets (Yaari et al., Frontiers in Immunology, 2013). "wrc" allows for only the WRC mutational hotspots to be included (where W equals A or T and R equals A or G). Specifying "all" will use all four types of mutations during SHM branching events, where the weights for each can be specified in the "SHM.nuc.prob" parameter.
#' @param SHM.nuc.prob Specifies the rate at which nucleotides change during speciation (SHM) events. This parameter depends on the type of mutation specified by SHM.method. For both "poisson" and "data", the input value determines the probability for each site to mutate (the whole sequence for "poisson" and the CDRs for "data"). For either "motif" or "wrc", the number of mutations per speciation event should be specified. Note that these are not probabilities, but the number of mutations that can occur (if the mutation is present in the sequence). If "all" is specified, the input should be a vector where the first element controls the poisson style mutations, second controls the "data", third controls the "motif" and fourth controls the "wrc".
#' @param SHM.isotype.dependent If TRUE, somatic hypermutation of certain isotype will happen based on probability specified in dataframe "iso_SHM_prob".
#' @param SHM.phenotype.dependent If TRUE, somatic hypermutation of certain phenotype will happen based on probability specified in dataframe "pheno_SHM_prob".
#' @param max.cell.number Integer value describing maximum number of cells allowed. Default is 1500.
#' @param max.clonotype.number Integer value describing maximum number of clones allowed. cell derived from the same mother cell belong to same clone.
#' @param death.rate Probability of cell death happen to each cell in each time step.
#' @param igraph.on If TRUE, mutational network for every B cell clone will be in the output. If False, the igraphs will not be included.
#' @param transcriptome.on If TRUE, the simulation will include transcriptome data. If FALSE, only vdj sequence will be simulated.
#' @param transcriptome.switch.independent TRUE or FALSE value describing whether transcriptome state is allowed to switch independently, not dependent on class switching or somatic hypermutation. If TURE, transcriptome.switch.prob should be specified to control the probability of transcriptome state switching.
#' @param transcriptome.switch.prob Probability of transcriptome state switching independently. Default values are in the attaching matrix "trans_switch_prob_b" and "trans_switch_prob_t". The order of cell type in rows and columns should be the same, and the order of the cell type in the matrix should match cell type names in transcriptome.states.
#' @param transcriptome.switch.isotype.dependent TRUE or FALSE value describing whether transcriptome state of a cell is allowed to switch depending on isotype switching. If TRUE, transcriptome state will switch once class switching happens.
#' @param transcriptome.switch.SHM.dependent TRUE or FALSE value describing whether transcriptome state of a cell is allowed to switch depending on somatic hypermutation. If TRUE, transcriptome state will switch once somatic hypermutation happens.
#' @param transcriptome.switch.selection.dependent If TRUE, selected cells will undergo transcriptome state switching if their transcriptome state is 1.
#' @param transcriptome.states A data.frame specifying base gene expression for different cell type, with gene names as row names, cell type names as column names. When missing, a default data.frame will be used. Default data.frame includes "GerminalcenterBcell", "NaiveBcell", "Plasmacell", "MemoryBcell" for B cells, and "Naïve Cd4", "ActivatedCd4", "MemoryCd4", "NaiveCd8", "EffectorCd8" ,"MemoryCd8","ExhaustedCd8" for T cells. The order of the cell type names in transcriptome.states should match cell type names in the transcriptome.switch.prob matrix.
#' @param transcriptome.noise A character expression specifying the distribution of noise ratio to be multiplied with the base gene expression for each cell. It should be a text expression that generates a numeric vector, which is of the same length as gene names in the trasncriptome.state input. Default value is "rnorm(nrow(transcriptome.states), mean = 1, sd = 0.3)".
#' @param seq.name Integer specifies how many top-ranking clones are included in Seq_Name dataframe in the output list for phylogenetic tree plotting in other pipeline. If missing, Seq_Name won't be included in the output.
#' @return A list containing the VDJ sequence and corresponding transcriptome data: "all_contig_annotations", "clonotypes", "all_contig", "consensus","reference","reference_real", "transcriptome","igraph_list_iso","igraph_list_trans","Seq_Name","igraph.index.attr","history","igraph.index","selected.seq","version","parameters".
#'
#' @export
simulate_repertoire <- function(initial.size.of.repertoire,
                                species,
                                cell.type,
                                cd4.proportion,
                                duration.of.evolution,
                                complete.duration,
                                vdj.productive,
                                vdj.model,
                                vdj.insertion.mean,
                                vdj.insertion.stdv,
                                vdj.branch.prob,
                                clonal.selection,
                                cell.division.prob,
                                sequence.selection.prob,
                                special.v.gene,
                                class.switch.prob,
                                class.switch.selection.dependent,
                                class.switch.independent,
                                SHM.method,
                                SHM.nuc.prob,
                                SHM.isotype.dependent,
                                SHM.phenotype.dependent,
                                max.cell.number,
                                max.clonotype.number,
                                death.rate,
                                igraph.on,
                                transcriptome.on,
                                transcriptome.switch.independent,
                                transcriptome.switch.prob,
                                transcriptome.switch.isotype.dependent,
                                transcriptome.switch.SHM.dependent,
                                transcriptome.switch.selection.dependent,
                                transcriptome.states,
                                transcriptome.noise,
                                seq.name
){

  # load polybox

  load(url("https://polybox.ethz.ch/index.php/s/zETU3ruyfTjj8T8/download"))

  # for cran check, global bind data sets----

  class_switch_prob_hum <- class_switch_prob_hum
  class_switch_prob_mus <- class_switch_prob_mus
  hum_b_h <- hum_b_h
  hum_b_l <- hum_b_l
  hum_b_trans <- hum_b_trans
  hum_t_h <- hum_t_h
  hum_t_l <- hum_t_l
  hum_t_trans <- hum_t_trans
  iso_SHM_prob <- iso_SHM_prob
  mus_b_h <- mus_b_h
  mus_b_l <- mus_b_l
  mus_b_trans <- mus_b_trans
  mus_t_h <- mus_t_h
  mus_t_l <- mus_t_l
  mus_t_trans <- mus_t_trans
  pheno_SHM_prob <- pheno_SHM_prob
  productive_seq <- productive_seq
  special_v <- special_v
  trans_switch_prob_b <- trans_switch_prob_b
  trans_switch_prob_t <- trans_switch_prob_t
  vae_seq <- vae_seq


  if(missing(species)) species<-"mus"
  if(missing(cell.type)) cell.type<-"B"
  if(missing(cd4.proportion)) cd4.proportion<-1
  if(missing(vdj.productive)) vdj.productive<-"random"
  if(missing(initial.size.of.repertoire)) initial.size.of.repertoire <- 10
  if(missing(vdj.model)) vdj.model <- "naive"
  if(missing(vdj.insertion.mean))  vdj.insertion.mean<-0.1
  if(missing(vdj.insertion.stdv))  vdj.insertion.stdv<-0.05
  if(missing(duration.of.evolution))  duration.of.evolution<-20
  if(missing(complete.duration)) complete.duration<-T
  if(missing(vdj.branch.prob))  vdj.branch.prob<-0.2
  if(missing(clonal.selection)) clonal.selection<-F
  if(missing(cell.division.prob))  cell.division.prob<-c(0.1,0.1)
  if(missing(sequence.selection.prob)) sequence.selection.prob<-0.01
  if(missing(special.v.gene)) special.v.gene<-F
  if(missing(class.switch.prob)&species=="mus")  class.switch.prob<-class_switch_prob_mus
  if(missing(class.switch.prob)&species=="hum")  class.switch.prob<-class_switch_prob_hum
  if(missing(class.switch.selection.dependent)) class.switch.selection.dependent<-F
  if(missing(class.switch.independent)) class.switch.independent<-T
  if(missing(SHM.nuc.prob))  SHM.nuc.prob<-0.001
  if(missing(SHM.method))  SHM.method<-"naive"
  if(missing(SHM.isotype.dependent)) SHM.isotype.dependent<-F
  if(missing(SHM.phenotype.dependent)) SHM.phenotype.dependent<-F
  if(missing(max.cell.number))  max.cell.number<-1500
  if(missing(max.clonotype.number))  max.clonotype.number<-20
  if(missing(death.rate)) death.rate<-0.001
  if(missing(igraph.on)&cell.type=="T") igraph.on<-F
  if(missing(igraph.on)&cell.type=="B") igraph.on<-T
  if(missing(transcriptome.on)) transcriptome.on<-T
  if(missing(transcriptome.switch.independent)) transcriptome.switch.independent<-T
  if(missing(transcriptome.switch.prob)&cell.type=="B") transcriptome.switch.prob<-trans_switch_prob_b
  if(missing(transcriptome.switch.prob)&cell.type=="T") transcriptome.switch.prob<-trans_switch_prob_t
  if(missing(transcriptome.switch.isotype.dependent)) transcriptome.switch.isotype.dependent<-F
  if(missing(transcriptome.switch.SHM.dependent)) transcriptome.switch.SHM.dependent<-F
  if(missing(transcriptome.switch.selection.dependent)) transcriptome.switch.selection.dependent<-F
  if(missing(transcriptome.states)&species=="mus"&cell.type=="B") transcriptome.states<-mus_b_trans
  if(missing(transcriptome.states)&species=="hum"&cell.type=="B") transcriptome.states<-hum_b_trans
  if(missing(transcriptome.states)&species=="mus"&cell.type=="T") transcriptome.states<-mus_t_trans
  if(missing(transcriptome.states)&species=="hum"&cell.type=="T") transcriptome.states<-hum_t_trans
  if(missing(transcriptome.noise)) transcriptome.noise<-"rnorm(nrow(transcriptome.states), mean = 1, sd = 0.3)"
  if(missing(seq.name)) seq.name<-F

  if(is.logical(transcriptome.switch.independent) & is.logical(transcriptome.switch.isotype.dependent) &is.logical( transcriptome.switch.SHM.dependent)==F){
    stop("transcriptome.switch.* is a T or F input")
  }
  if(!(nrow(transcriptome.states)==length(eval(parse(text=transcriptome.noise))))){
    stop("transcriptome.noise vector length different from transcriptome gene length")
  }
  if(length(transcriptome.noise)>1&length(transcriptome.noise)!=ncol(transcriptome.states)){
    stop("plase assgin as many transcriptome.noise vector items as the number of transcriptome.states or just assign a common one")
  }
  if(!missing(transcriptome.states)) cd4.proportion<-1

#check the order of cell names in transcriptome.switch.prob and transcriptome.states
 if(all(!(rownames(transcriptome.switch.prob) == colnames(transcriptome.switch.prob)&rownames(transcriptome.switch.prob) == colnames(transcriptome.states)))){
   stop("transcriptome.state colnames and transcriptome.switch.prob colnames and rownames should be the same")
 }

 if(length(cell.division.prob)==1){
   cell.division.prob<-c(cell.division.prob,cell.division.prob)
 }
  if(max(cell.division.prob)>1){
    stop("Any value in cell.division.prob should less than 1")
  }


# require( stringr, reshape2, dplyr, do, igraph, stats, ggplot2)

#initial vectors------
  barcode<-c()
  raw_contig_id<-c()
  chain<-c()
  v_gene<-c()
  d_gene<-c()
  j_gene<-c()
  c_gene<-c()
  v_seq<-c()
  d_seq<-c()
  j_seq<-c()

  ref_seq<-c()
  raw_clonotype_id<-c()
  clonotype_num<-initial.size.of.repertoire

  seq<-c()
  Length<-c()
  #single format
  clonotype_id<-c()
  seq_combi<-c()
  barcode_uniq<-c()
  gen<-c()
  isotype<-c()
  seq.number<-c()
  seq.history<-c()
  igraph.index<-list()
  barcode.history<-c()
  igraph.index.attr<-list()
  igraph_list<-list()
  igraph.index.jr<-list()
  colors<-colors
  pie.values.list<-list()
  pie.values<-list()
  trans_dis<-c()
  trans_ls<-list()
  trans_state<-c()
  trans_state_his<-c()

  igraph_list_trans<-list()
  pie_values_trans_list<-list()
  cd4_prob<-cd4.proportion

  cdr3<-c()
  cdr3_nt<-c()

  selected_seq<-c()
  selected_rate<-c()

#load input----
  if(length(transcriptome.noise)>1){

    for (i in 1:ncol(transcriptome.states)){
      trans_dis[i]<-paste0("transcriptome.states[",i,"]*",transcriptome.noise[i])
    }
  }
  if(length(transcriptome.noise)==1){
    for (i in 1:ncol(transcriptome.states)){
      trans_dis[i]<-paste0("transcriptome.states[",i,"]*",transcriptome.noise)
    }
  }

  if(cell.type=="B" & species=="mus" | species=="mouse" | species=="blc6"){
  if(vdj.productive=="random"){
    seq_input_h<-mus_b_h
    seq_input_l<-mus_b_l
    }
  if(vdj.productive=="naive"){
    seq_input_h<-productive_seq$mus_b_h
    seq_input_l<-productive_seq$mus_b_l
    # ref_seq_h<-mus_b_h
    # ref_seq_l<-mus_b_l
  }
  if(vdj.productive=="vae"){
    seq_input_h<-vae_seq$mus_b_h
    seq_input_l<-vae_seq$mus_b_l
    # ref_seq_h<-mus_b_h
    # ref_seq_l<-mus_b_l
  }
    c_genename_h<-"IGHM"
  }
  if(cell.type=="B" & species=="hum" | species=="human"){
    if(vdj.productive=="random"){
    seq_input_h<-hum_b_h
    seq_input_l<-hum_b_l
    }
    if(vdj.productive=="naive"){
    seq_input_h<-productive_seq$hum_b_h
    seq_input_l<-productive_seq$hum_b_l
    # ref_seq_h<-hum_b_h
    # ref_seq_l<-hum_b_l
    }
    c_genename_h<-"IGHM"
  if(vdj.productive=="vae"){
    seq_input_h<-vae_seq$hum_b_h
    seq_input_l<-vae_seq$hum_b_l
    # ref_seq_h<-hum_b_h
    # ref_seq_l<-hum_b_l
  }
  c_genename_h<-"IGHM"
}
  if(cell.type=="T" & species=="mus" | species=="mouse" | species=="blc6" ){
    if(vdj.productive=="random"){
    seq_input_h<-mus_t_h
    seq_input_l<-mus_t_l
    }
    if(vdj.productive=="naive"){
    seq_input_h<-productive_seq$mus_t_h
    seq_input_l<-productive_seq$mus_t_l
    # ref_seq_h<-mus_t_h
    # ref_seq_l<-mus_t_l
    }
    if(vdj.productive=="vae"){
      seq_input_h<-vae_seq$mus_t_h
      seq_input_l<-vae_seq$mus_t_l
      # ref_seq_h<-mus_t_h
      # ref_seq_l<-mus_t_l
    }
    c_genename_h<-"TRBC"
  }
  if(cell.type=="T" & species=="hum" | species=="human"){
    if(vdj.productive=="random"){
    seq_input_h<-hum_t_h
    seq_input_l<-hum_t_l
    }
    if(vdj.productive=="naive"){
    seq_input_h<-productive_seq$hum_t_h
    seq_input_l<-productive_seq$hum_t_l
    # ref_seq_h<-hum_t_h
    # ref_seq_l<-hum_t_l
    }
    if(vdj.productive=="vae"){
      seq_input_h<-vae_seq$hum_t_h
      seq_input_l<-vae_seq$hum_t_l
      # ref_seq_h<-hum_t_h
      # ref_seq_l<-hum_t_l
    }
    c_genename_h<-"TRBC"
  }



#generate initial barcode vector------
barcode<-.GENERATE.BARCODE(initial.size.of.repertoire)

#initial----
for(i in 1:initial.size.of.repertoire){
  if(vdj.productive=="random"){
  #heavy chain seq
  indV <- sample(x = 1:nrow(seq_input_h[[1]]),size = 1,replace = FALSE)
  indD <- sample(x = 1:nrow(seq_input_h[[2]]),size = 1,replace=FALSE)
  indJ <- sample(x = 1:nrow(seq_input_h[[3]]),size = 1,replace=FALSE)
  #light chain seq
  indVL <- sample(x = 1:nrow(seq_input_l[[1]]),size = 1,replace = FALSE)
  indJL <- sample(x = 1:nrow(seq_input_l[[2]]),size=1,replace=FALSE)
  #seq
  seq[2*i-1]<- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(seq_input_h[[1]][[2]][indV]),as.character(seq_input_h[[2]][[2]][indD]),as.character(seq_input_h[[3]][[2]][indJ]),
                                                   method=vdj.model,
                                                   chain.type="heavy",
                                                   species=species,
                                                   vdj.insertion.mean=vdj.insertion.mean,
                                                   vdj.insertion.stdv=vdj.insertion.stdv))
  seq[2*i] <- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(seq_input_l[[1]][[2]][indVL]),"",as.character(seq_input_l[[2]][[2]][indJL]),
                                                  method=vdj.model,
                                                  chain.type="light",
                                                  species=species,
                                                  vdj.insertion.mean=vdj.insertion.mean,
                                                  vdj.insertion.stdv=vdj.insertion.stdv))

  #vdjc gene names
  v_gene[2*i-1]<-as.character(seq_input_h[[1]][[1]][indV])
  v_gene[2*i]<-as.character(seq_input_l[[1]][[1]][indVL])
  d_gene[2*i-1]<-as.character(seq_input_h[[2]][[1]][indD])
  d_gene[2*i]<-"None"
  j_gene[2*i-1]<-as.character(seq_input_h[[3]][[1]][indJ])
  j_gene[2*i]<-as.character(seq_input_l[[2]][[1]][indJL])
  c_gene[2*i-1]<-c_genename_h
  c_gene[2*i]<-paste0(stringr::str_sub(v_gene[2*i],1,3),"C")
  #ref seq
  v_seq[2*i-1]<-as.character(seq_input_h[[1]][[2]][indV])
  v_seq[2*i]<-as.character(seq_input_l[[1]][[2]][indVL])
  d_seq[2*i-1]<-as.character(seq_input_h[[2]][[2]][indD])
  d_seq[2*i]<-"None"
  j_seq[2*i-1]<-as.character(seq_input_h[[3]][[2]][indJ])
  j_seq[2*i]<-as.character(seq_input_l[[2]][[2]][indJL])
  #the very first ref seq
  ref_seq[2*i-1]<-seq[2*i-1]
  ref_seq[2*i]<-seq[2*i]
  }
  if(vdj.productive=="naive"|vdj.productive=="vae"){
    indH<-sample(x=1:nrow(seq_input_h),size=1)
    indL<-sample(x=1:nrow(seq_input_l),size=1)
    #seq
    seq[2*i-1]<- tolower(seq_input_h$targetSequences[indH])
    seq[2*i] <- tolower(seq_input_l$targetSequences[indL])
    #ref seq
    v_seq[2*i-1]<-tolower(seq_input_h$v_seq[indH])
    v_seq[2*i]<-tolower(seq_input_l$v_seq[indL])
    #which
    indH<-which(toupper(seq[2*i-1])==seq_input_h$targetSequences)
    indL<-which(toupper(seq[2*i])==seq_input_l$targetSequences)
    #vdjc gene names
    v_gene[2*i-1]<-seq_input_h$bestVHit[indH]
    v_gene[2*i]<-seq_input_l$bestVHit[indL]
    d_gene[2*i-1]<-seq_input_h$bestDHit[indH]
    d_gene[2*i]<-"None"
    j_gene[2*i-1]<-seq_input_h$bestJHit[indH]
    j_gene[2*i]<-seq_input_l$bestJHit[indL]
    c_gene[2*i-1]<-c_genename_h
    c_gene[2*i]<-paste0(stringr::str_sub(v_gene[2*i],1,3),"C")

    #cdr3
    cdr3[2*i-1]<-seq_input_h$aaSeqCDR3[indH]
    cdr3[2*i]<-seq_input_h$aaSeqCDR3[indL]
    cdr3_nt[2*i-1]<-seq_input_h$nSeqCDR3[indH]
    cdr3_nt[2*i]<-seq_input_h$nSeqCDR3[indL]
    #ref seq
    ref_seq[2*i-1]<-seq[2*i-1]
    ref_seq[2*i]<-seq[2*i]
  }

  #chain
  chain[2*i-1]<-stringr::str_sub(v_gene[2*i-1],1,3)
  chain[2*i]<-stringr::str_sub(v_gene[2*i],1,3)
  #length
  Length[2*i-1]<-nchar(seq[2*i-1])
  Length[2*i]<-nchar(seq[2*i])
  #contig
  raw_contig_id[2*i-1]<-paste(barcode[2*i-1],"contig_1",sep="_")
  raw_contig_id[2*i]<-paste(barcode[2*i],"contig_2",sep="_")
  #clonotype
  raw_clonotype_id[2*i-1]<-paste0("clonotype",i)
  raw_clonotype_id[2*i]<-paste0("clonotype",i)
  #single
  seq_combi[i]<- paste0(chain[2*i-1],":",seq[2*i-1],";",chain[2*i],":",seq[2*i])
  clonotype_id[i]<- paste0("clonotype",i)
  barcode_uniq[i]<-barcode[2*i]
  gen[i]<-1

  seq.number[i]<-i
  seq.history[i]<-seq_combi[i]
  igraph.index[[i]]<-c(seq.number[i])
  isotype[i]<-1

  trans_state[i]<-1
  if(cell.type=="T"&length(cd4_prob)>0){
    trans_state[i]<-sample(x=c(1,4),size = 1,replace = T, prob = c(cd4_prob,1-cd4_prob))
  }
  trans_state_his[i]<- trans_state[i]

  #special v_h&v_l combi check----

  if(length(cell.division.prob)>2&special.v.gene==T){

    pei<-which(special_v[,1]==v_gene[2*i-1])
    pei2<-which(special_v[,2]==v_gene[2*i])
    exe<-F
    if(length(pei)>0){
      if(length(pei2)>0|any(special_v[pei,2]=="")) exe<-T
    }
    if(length(pei2)>0){
      if(length(pei)>0|any(special_v[pei2,1]=="")) exe<-T
    }
    if(exe==T){
      special_row<-c(pei,pei2)[which(c(pei,pei2)>0)[1]]
      seq_selection<-sample(x=c(0,1),size = 1,replace = T, prob = c(special_v[,3][special_row],1-special_v[,3][special_row]))
      if(seq_selection==0){
        selected_rate<-c(selected_rate,stats::runif(1,cell.division.prob[3],cell.division.prob[length(cell.division.prob)]))
        selected_seq<-c(selected_seq,seq.number[i])
      }
    }
  }
}

  #check sequence selection initial----
  if(length(cell.division.prob)>2){
    seq_selection<-sample(x=c(0,1),size = length(seq.number),replace = T, prob = c(sequence.selection.prob,1-sequence.selection.prob))
    selected_rate<-c(selected_rate,stats::runif(sum(seq_selection==0),cell.division.prob[3],cell.division.prob[length(cell.division.prob)]))
    selected_seq<-c(selected_seq,seq.number[which(seq_selection==0)])
    if(special.v.gene==T){
    selected_rate<-selected_rate[!duplicated(selected_seq)]
    selected_seq<-selected_seq[!duplicated(selected_seq)]
    }
  }

  barcode.history<-barcode_uniq
#evolution----
  if(duration.of.evolution>=1){
    for(i in 1:duration.of.evolution){
  #new vdj----
    is_new_VDJ <- sample(x=c(0,1), replace=TRUE,size = 1, prob=c(vdj.branch.prob,1- vdj.branch.prob))
    if (is_new_VDJ==0 & length(barcode_uniq)<max.cell.number & clonotype_num< max.clonotype.number){
      clonotype_num<-clonotype_num+1
      l<-length(barcode)
      barcode<-.ADD.BARCODE(barcode)

      if(vdj.productive=="random"){
      #heavy chain seq
      indV <- sample(x = 1:nrow(seq_input_h[[1]]),size = 1,replace = FALSE)
      indD <- sample(x = 1:nrow(seq_input_h[[2]]),size = 1,replace=FALSE)
      indJ <- sample(x = 1:nrow(seq_input_h[[3]]),size = 1,replace=FALSE)
      #light chain seq
      indVL <- sample(x = 1:nrow(seq_input_l[[1]]),size = 1,replace = FALSE)
      indJL <- sample(x = 1:nrow(seq_input_l[[2]]),size=1,replace=FALSE)
      #seq
      seq[l+1]<- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(seq_input_h[[1]][[2]][indV]),as.character(seq_input_h[[2]][[2]][indD]),as.character(seq_input_h[[3]][[2]][indJ]),
                                                     method=vdj.model,
                                                     chain.type="heavy",
                                                     species=species,
                                                     vdj.insertion.mean=vdj.insertion.mean,
                                                     vdj.insertion.stdv=vdj.insertion.stdv))
      seq[l+2] <- as.character(.VDJ_RECOMBIN_FUNCTION(as.character(seq_input_l[[1]][[2]][indVL]),"",as.character(seq_input_l[[2]][[2]][indJL]),
                                                      method=vdj.model,
                                                      chain.type="light",
                                                      species=species,
                                                      vdj.insertion.mean=vdj.insertion.mean,
                                                      vdj.insertion.stdv=vdj.insertion.stdv))



      #vdjc gene names
      v_gene[l+1]<-as.character(seq_input_h[[1]][[1]][indV])
      v_gene[l+2]<-as.character(seq_input_l[[1]][[1]][indVL])
      d_gene[l+1]<-as.character(seq_input_h[[2]][[1]][indD])
      d_gene[l+2]<-"None"
      j_gene[l+1]<-as.character(seq_input_h[[3]][[1]][indJ])
      j_gene[l+2]<-as.character(seq_input_l[[2]][[1]][indJL])
      c_gene[l+1]<-c_genename_h
      c_gene[l+2]<-paste0(stringr::str_sub(v_gene[l+2],1,3),"C")

      #vdj ref seq
      v_seq[l+1]<-as.character(seq_input_h[[1]][[2]][indV])
      v_seq[l+2]<-as.character(seq_input_l[[1]][[2]][indVL])
      d_seq[l+1]<-as.character(seq_input_h[[2]][[2]][indD])
      d_seq[l+2]<-"None"
      j_seq[l+1]<-as.character(seq_input_h[[3]][[2]][indJ])
      j_seq[l+2]<-as.character(seq_input_l[[2]][[2]][indJL])
      #another ref seq
      r<-length(ref_seq)
      ref_seq[r+1]<-seq[l+1]
      ref_seq[r+2]<-seq[l+2]
      }
      if(vdj.productive=="naive"|vdj.productive=="vae"){
        indH<-sample(x=1:nrow(seq_input_h),size=1)
        indL<-sample(x=1:nrow(seq_input_l),size=1)
        #seq
        seq[l+1]<- tolower(seq_input_h$targetSequences[indH])
        seq[l+2] <- tolower(seq_input_l$targetSequences[indL])
        #ref seq
        v_seq[l+1]<-tolower(seq_input_h$v_seq[indH])
        v_seq[l+2]<-tolower(seq_input_l$v_seq[indL])
        #vdjc gene names
        v_gene[l+1]<-seq_input_h$bestVHit[indH]
        v_gene[l+2]<-seq_input_l$bestVHit[indL]
        d_gene[l+1]<-seq_input_h$bestDHit[indH]
        d_gene[l+2]<-"None"
        j_gene[l+1]<-seq_input_h$bestJHit[indH]
        j_gene[l+2]<-seq_input_l$bestJHit[indL]
        c_gene[l+1]<-c_genename_h
        c_gene[l+2]<-paste0(stringr::str_sub(v_gene[l+2],1,3),"C")
        #cdr3
        cdr3[l+1]<-seq_input_h$aaSeqCDR3[indH]
        cdr3[l+2]<-seq_input_l$aaSeqCDR3[indL]
        cdr3_nt[l+1]<-seq_input_h$nSeqCDR3[indH]
        cdr3_nt[l+2]<-seq_input_l$nSeqCDR3[indL]
        #ref seq
        r<-length(ref_seq)
        ref_seq[r+1]<-seq[l+1]
        ref_seq[r+2]<-seq[l+2]
     }
    #common updates
      #chain
      chain[l+1]<-stringr::str_sub(v_gene[l+1],1,3)
      chain[l+2]<-stringr::str_sub(v_gene[l+2],1,3)
      #length
      Length[l+1]<-nchar(seq[l+1])
      Length[l+2]<-nchar(seq[l+2])
      #contig
      raw_contig_id[l+1]<-paste(barcode[l+1],"contig_1",sep="_")
      raw_contig_id[l+2]<-paste(barcode[l+2],"contig_2",sep="_")
      #clonotype
      raw_clonotype_id[l+1]<-paste("clonotype",clonotype_num,sep="")
      raw_clonotype_id[l+2]<-paste("clonotype",clonotype_num,sep="")
      #single format used for finding consensus later
      seq_combi[0.5*l+1]<- paste(chain[l+1],":",seq[l+1],";",chain[l+2],":",seq[l+2],sep="")
      clonotype_id[0.5*l+1]<- raw_clonotype_id[l+1]
      barcode_uniq[0.5*l+1]<- barcode[l+2]
      gen[0.5*l+1]<-1
      trans_state[0.5*l+1]<-1
      #history update
      hisl<-length(seq.history)
      barcode.history[hisl+1]<-barcode_uniq[0.5*l+1]
      seq.number[hisl+1]<-max(seq.number)+1
      seq.history[hisl+1]<-seq_combi[0.5*l+1]
      igraph.index[[clonotype_num]]<-c(seq.number[hisl+1])
      isotype[hisl+1]<-1
      ####!!!!
      if(cell.type=="T"&length(cd4_prob)>0){
        trans_state[0.5*l+1]<-sample(x=c(1,4),size = 1,replace = T, prob = c(cd4_prob,1-cd4_prob))
      }
      trans_state_his[hisl+1]<-trans_state[0.5*l+1]
      #check sequence selection in new VDJ####
      if(length(cell.division.prob)>2){
        seq_selection<-sample(x=c(0,1),size = 1, prob = c(sequence.selection.prob,1-sequence.selection.prob))

        if(special.v.gene==T){
          pei<-which(special_v[,1]==v_gene[l+1])
          pei2<-which(special_v[,2]==v_gene[l+2])
          exe<-F
          if(length(pei)>0){
            if(length(pei2)>0|any(special_v[pei,2]=="")) exe<-T
          }
          if(length(pei2)>0){
            if(length(pei)>0|any(special_v[pei2,1]=="")) exe<-T
          }
          if(exe==T){
            special_row<-c(pei,pei2)[which(c(pei,pei2)>0)[1]]
            seq_selection<-sample(x=c(0,1),size = 1,replace = T, prob = c(special_v[,3][special_row],1-special_v[,3][special_row]))
            if(seq_selection==0){
              selected_rate<-c(selected_rate,stats::runif(1,cell.division.prob[3],cell.division.prob[length(cell.division.prob)]))
              selected_seq<-c(selected_seq,seq.number[hisl+1])
            }
          }
        }
        if(special.v.gene==F){
          if(seq_selection==0){
            selected_rate<-c(selected_rate,stats::runif(1,cell.division.prob[3],cell.division.prob[length(cell.division.prob)]))
            selected_seq<-c(selected_seq,seq.number[hisl+1])
          }
        }
      }

    }
    if(complete.duration==F&!(length(barcode_uniq)<max.cell.number & clonotype_num< max.clonotype.number)) break
  #cell division----
    if(length(barcode_uniq)<max.cell.number & clonotype_num< max.clonotype.number){
      #check existing cell division for every cell
      for(j in 1:length(barcode_uniq)){

        if(clonal.selection==F){
        is_cell_division<- .CELL.DIVISION.LINEAR.INVERSE(clonotype_id,j,cell.division.prob)
        }

        if(clonal.selection==T){
        is_cell_division<- .CELL.DIVISION.LINEAR(clonotype_id,j,cell.division.prob)
        }

        if(length(cell.division.prob)>2){
          seq_number_j<-seq.number[which(barcode.history==barcode_uniq[j])]
          if(length(which(selected_seq==seq_number_j))>0){
            rate_j<-selected_rate[which(selected_seq==seq_number_j)]
            #!!!!!wenti
            is_cell_division<-sample(x=c(0,1),size=1,prob =c(rate_j,1-rate_j))
          }
        }
        if (is_cell_division==0&length(barcode_uniq)<max.cell.number & clonotype_num< max.clonotype.number ){

          l<-length(barcode)
          #barcode
          barcode<-.ADD.BARCODE(barcode)
          #seq stays the same
          seq[l+1]<-seq[2*j-1]
          seq[l+2]<-seq[2*j]
          #length
          Length[l+1]<-nchar(seq[2*j-1])
          Length[l+2]<-nchar(seq[2*j])
          #contig
          raw_contig_id[l+1]<-paste(barcode[l+1],"contig_1",sep="_")
          raw_contig_id[l+2]<-paste(barcode[l+2],"contig_2",sep="_")
          #clonotype stay the same
          raw_clonotype_id[l+1]<-clonotype_id[j]
          raw_clonotype_id[l+2]<-clonotype_id[j]
          #vdj seq for reference
          if(vdj.productive=="random"){
            d_seq[l+1]<-d_seq[2*j-1]
            d_seq[l+2]<-"None"
            j_seq[l+1]<-j_seq[2*j-1]
            j_seq[l+2]<-j_seq[2*j]
          }

          v_seq[l+1]<-v_seq[2*j-1]
          v_seq[l+2]<-v_seq[2*j]
          #gene names
          v_gene[l+1]<-v_gene[2*j-1]
          v_gene[l+2]<-v_gene[2*j]
          d_gene[l+1]<-d_gene[2*j-1]
          d_gene[l+2]<-"None"
          j_gene[l+1]<-j_gene[2*j-1]
          j_gene[l+2]<-j_gene[2*j]
          c_gene[l+1]<-c_gene[2*j-1]
          c_gene[l+2]<-c_gene[2*j]
          #cdr3
          if(vdj.productive=="naive"|vdj.productive=="vae"){
            cdr3[l+1]<-cdr3[2*j-1]
            cdr3[l+2]<-cdr3[2*j]
            cdr3_nt[l+1]<-cdr3_nt[2*j-1]
            cdr3_nt[l+2]<-cdr3_nt[2*j]
          }
          #chain names
          chain[l+1]<-stringr::str_sub(v_gene[l+1],1,3)
          chain[l+2]<-stringr::str_sub(v_gene[l+2],1,3)
          #single format
          trans_state[0.5*l+1]<-trans_state[j]
          seq_combi[0.5*l+1]<- paste(chain[l+1],":", seq[l+1],";",chain[l+2],":",seq[l+2],sep ="")
          clonotype_id[0.5*l+1]<- raw_clonotype_id[l+1]
          barcode_uniq[0.5*l+1]<-barcode[l+1]
          gen[0.5*l+1]<-gen[j]+1

          #history format
          n<- which(barcode.history==barcode_uniq[j])
          hisl<-length(seq.history)

          isotype[hisl+1]<-isotype[n]
          barcode.history[hisl+1]<-barcode[l+1]
          seq.number[hisl+1]<-seq.number[n]#seq number should be the same as the mother cell
          seq.history[hisl+1]<-seq.history[n]
          trans_state_his[hisl+1]<-trans_state_his[n]

        }
      }# for loop end
    }#if number check end, cell division end
    if(complete.duration==F&!(length(barcode_uniq)<max.cell.number & clonotype_num< max.clonotype.number)) break
  #transcriptome switch independent----
    trans_switch<-rep(F,length(barcode_uniq))
    if(transcriptome.switch.independent== T){
      for(k in 1: length(barcode_uniq)) {
        switched<-.TRANS.SWITCH(trans_state[k],transcriptome.switch.prob)
        if(switched[[2]]){
        trans_state[k]<-switched[[1]]
        m<-which(barcode.history==barcode_uniq[k])
        trans_state_his[m]<-switched[[1]]
        trans_switch[k]<-T
        }
      }
    }
    if(transcriptome.switch.selection.dependent==T){
      for(k in 1: length(barcode_uniq)) {
        m<-which(barcode.history==barcode_uniq[k])
        if(length(which(selected_seq==seq.number[m]))>0&trans_state[k]==1){
          switched<-.TRANS.SWITCH.DEPENDANT(trans_state[k],transcriptome.switch.prob)
          trans_state[k]<-switched[[1]]
          m<-which(barcode.history==barcode_uniq[k])
          trans_state_his[m]<-switched[[1]]
          trans_switch[k]<-T
        }
      }
    }

  #b cell special

  #class switch & mutation----
    if (cell.type =="B"){
      for(k in 1: length(barcode_uniq)) {
        # class switch----
        n<- which(barcode.history==barcode_uniq[k])
        selected<-(length(which(selected_seq==seq.number[n]))>0)
        if(class.switch.independent==T){
        switched<-.TRANS.SWITCH(isotype[n],class.switch.prob)
          if(switched[[2]]){
          isotype[n]<-switched[[1]]
          c_gene[2*k-1]<-colnames(class.switch.prob)[switched[[1]]]
          }
        }
        #selected sequence make sure switch
        if(class.switch.selection.dependent==T&isotype[n]==1&selected){
          switched<-.TRANS.SWITCH.DEPENDANT(isotype[n],class.switch.prob)
          isotype[n]<-switched[[1]]
          c_gene[2*k-1]<-colnames(class.switch.prob)[switched[[1]]]
        }
        #trans switch isotype dependent----
        if(switched[[2]]&trans_switch[k] == F & transcriptome.switch.isotype.dependent== T&sum(transcriptome.switch.prob[trans_state[k],])>0){
            switched<-.TRANS.SWITCH.DEPENDANT(trans_state[k],transcriptome.switch.prob)
            trans_state[k]<-switched[[1]]
            m<-which(barcode.history==barcode_uniq[k])
            trans_state_his[m]<-switched[[1]]
            trans_switch[k]<-T
          }

        #mutation----

        if((SHM.isotype.dependent==F|length(which(iso_SHM_prob[,1]==c_gene[2*k-1]))==0)&
           (SHM.phenotype.dependent==F|length(which(pheno_SHM_prob[,1]==colnames(transcriptome.switch.prob)[trans_state[k]]))==0)
           ){
          mut_h<-.SHM_FUNCTION_SEQUENCE4(seq[2*k-1],SHM.method,v_seq[2*k-1],SHM.nuc.prob)
          mut_l<-.SHM_FUNCTION_SEQUENCE4(seq[2*k],SHM.method,v_seq[2*k],SHM.nuc.prob)
          Ctrl<-mut_h[[2]]+mut_l[[2]]
        }

        if(!(
          (SHM.isotype.dependent==F|length(which(iso_SHM_prob[,1]==c_gene[2*k-1]))==0)&
          (SHM.phenotype.dependent==F|length(which(pheno_SHM_prob[,1]==colnames(transcriptome.switch.prob)[trans_state[k]]))==0)
           )){
          new.SHM.nuc.prob<-max(iso_SHM_prob[which(iso_SHM_prob[,1]==c_gene[2*k-1]),2],pheno_SHM_prob[which(pheno_SHM_prob[,1]==colnames(transcriptome.switch.prob)[trans_state[k]]),2])
          mut_h<-.SHM_FUNCTION_SEQUENCE4(seq[2*k-1],SHM.method,v_seq[2*k-1],new.SHM.nuc.prob)
          mut_l<-.SHM_FUNCTION_SEQUENCE4(seq[2*k],SHM.method,v_seq[2*k],new.SHM.nuc.prob)
          Ctrl<-mut_h[[2]]+mut_l[[2]]
        }


        if(Ctrl>0){
          #update seq
          seq[2*k-1]<- mut_h[[1]]
          seq[2*k]<-mut_l[[1]]

          #update length
          Length[2*k-1]<-nchar(seq[2*k-1])
          Length[2*k]<-nchar(seq[2*k])

          #update single format
          seq_combi[k]<- paste(chain[2*k-1],":",seq[2*k-1],";",chain[2*k],":",seq[2*k],sep ="")

          #update igraph history
          hisl<-length(seq.history)
          m<-which(barcode.history == barcode_uniq[k])

          barcode.history[m]<-NA
          barcode.history[hisl+1]<-barcode_uniq[k]
          seq.number[hisl+1]<-max(seq.number)+1
          seq.history[hisl+1]<-seq_combi[k]
          isotype[hisl+1]<-isotype[m]
          trans_state_his[hisl+1] <- trans_state_his[m]

          clonotype.number<-as.numeric(stringr::str_sub(clonotype_id[k],10,-1))
          igraph.index[[clonotype.number]]<-c(igraph.index[[clonotype.number]],seq.number[m],seq.number[hisl+1])

          #check sequence selection in mutation
          if(length(cell.division.prob)>2){
            seq_selection<-sample(x=c(0,1),size = 1,prob = c(sequence.selection.prob,1-sequence.selection.prob))
            if(seq_selection==0){
              selected_rate<-c(selected_rate,stats::runif(1,cell.division.prob[3],cell.division.prob[length(cell.division.prob)]))
              selected_seq<-c(selected_seq,seq.number[hisl+1])
            }
          }

        #try: modify cdr3
        #trans switch SHM dependent
        if(trans_switch[k] == F & transcriptome.switch.isotype.dependent== T&sum(transcriptome.switch.prob[trans_state[k],])>0){
          switched<-.TRANS.SWITCH.DEPENDANT(trans_state[k],transcriptome.switch.prob)
          trans_state[k]<-switched[[1]]
          m<-which(barcode.history==barcode_uniq[k])
          trans_state_his[m]<-switched[[1]]
          trans_switch[k]<-T
        }
      }
      }
    }#B cell special end
  #cell death----
    if(length(barcode_uniq<max.cell.number)&clonotype_num<max.clonotype.number){
      for(k in 1:length(barcode_uniq)){
        is_death <- sample(x=c(0,1), replace=TRUE,size = 1, prob=c(death.rate, 1- death.rate))
        if (is_death ==0){
          #fix row for double and history form first
          del<-c(-(2*k-1),-2*k)
          del.his<-which(barcode.history==barcode_uniq[k])
          #delete double format
          barcode<-barcode[del]
          seq<-seq[del]
          Length<-Length[del]
          raw_contig_id<-raw_contig_id[del]
          raw_clonotype_id<-raw_clonotype_id[del]
          v_gene<-v_gene[del]
          d_gene<-d_gene[del]
          j_gene<-j_gene[del]
          c_gene<-c_gene[del]
          chain<-chain[del]
          v_seq<- v_seq[del]
          d_seq<-d_seq[del]
          j_seq<-j_seq[del]
          cdr3<-cdr3[del]
          cdr3_nt<-cdr3_nt[del]
          #single format
          barcode_uniq<-barcode_uniq[-k]
          seq_combi<-seq_combi[-k]
          clonotype_id<-clonotype_id[-k]
          gen<-gen[-k]
          trans_state<-trans_state[-k]

          #history format
          barcode.history[del.his]<-NA
        }
      }
    }#cell death end
  }#evolution ends----
  }
#write data frames----
  #all_contig_annotations

  #dummy
  is_cell<-rep("Ture",length(barcode))
  high_confidence<-is_cell
  full_length<-is_cell
  productive<-is_cell
  reads<-rep(1000,length(barcode))
  umis<-rep(10,length(barcode))

  if(vdj.productive=="random"){
  cdr3<-rep("KKKK",length(barcode))
  cdr3_nt<-rep("AAAAAAAAAAAA",length(barcode))
  all_contig_annotations<-data.frame(barcode,is_cell,raw_contig_id,high_confidence,Length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,umis,raw_clonotype_id)
  all_contig_annotations$raw_consensus_id<-paste(raw_clonotype_id,"consensus",substr(raw_contig_id,nchar(raw_contig_id),nchar(raw_contig_id)),sep = "_")
  }

   #clonotypes from contig annotations
  clonotypes<-as.data.frame(table(clonotype_id))
  colnames(clonotypes)<-c("clonotype_id","frequency")
  clonotypes$proportion<-clonotypes$frequency/sum(clonotypes$frequency)

  if(vdj.productive=="naive"|vdj.productive=="vae"){
    #single format cdr3 for clonotype
    cdr3_combi<-.GENERATE.COMBI(cdr3,chain)
    cdr3_nt_combi<-.GENERATE.COMBI(cdr3_nt,chain)
    cdr3_df<-data.frame(clonotype_id,cdr3_combi,cdr3_nt_combi)
    clonotypes<-merge(clonotypes,cdr3_df,by="clonotype_id",all.x=T)
    clonotypes<-clonotypes[!(duplicated(clonotypes$clonotype_id)),]
    names(clonotypes)[4:5]<-c("cdr3","cdr3_nt")
    all_contig_annotations<-data.frame(barcode,is_cell,raw_contig_id,high_confidence,Length,chain,v_gene,d_gene,j_gene,c_gene,full_length,productive,cdr3,cdr3_nt,reads,raw_clonotype_id)
  }

  raw_consensus<-.WRITE.CONSENSUS(clonotype_id,sequence.combined = seq_combi,barcode.unique = barcode_uniq,clonotypes.dataframe = clonotypes)
  # all_contig_annotations$clonotype_chain<-paste(all_contig_annotations$raw_clonotype_id,stringr::str_sub(all_contig_annotations$raw_contig_id,-1,-1),sep="_")
  # all_contig_annotations<-merge(all_contig_annotations,raw_consensus,by="clonotype_chain",all=TRUE)
  # all_contig_annotations<-subset(all_contig_annotations,select=c("barcode","raw_contig_id","Length","chain","v_gene", "d_gene","j_gene","c_gene","raw_clonotype_id","raw_consensus_id"))
  consensus<-raw_consensus[,-1]


  #all_contig
  all_contig<-data.frame(raw_contig_id,seq)
  #reference
  if(vdj.productive=="random"){
    reference_id<-paste0(raw_clonotype_id,"_concat_ref_",stringr::str_sub(raw_contig_id,-1,-1))
    reference_seq<-c()
    d_seq1<-do::Replace(data=d_seq,from="None",to="")
    reference_seq<-paste0(v_seq,d_seq1,j_seq)
    reference<-data.frame(reference_id,reference_seq)
    reference<-reference[!(duplicated(reference$reference_id)),]

    reference_id_real<-c()
    for(cn in 1:clonotype_num){
      reference_id_real[2*cn-1]<-paste0("clonotype",cn,"_concat_ref_",1)
      reference_id_real[2*cn]<-paste0("clonotype",cn,"_concat_ref_",2)
    }
    reference_real<-data.frame(reference_id_real,ref_seq)
    colnames(reference_real)<-c("reference_id","reference_seq")
  }
  if(vdj.productive=="naive"|vdj.productive=="vae"){
    reference_id<-c()
    for(cn in 1:clonotype_num){
      reference_id[2*cn-1]<-paste0("clonotype",cn,"_concat_ref_",1)
      reference_id[2*cn]<-paste0("clonotype",cn,"_concat_ref_",2)
    }
    reference_real<-data.frame(reference_id,ref_seq)
    colnames(reference_real)<-c("reference_id","reference_seq")
    reference<-reference_real
  }
  history<-data.frame(seq.number,barcode.history,seq.history,isotype,trans_state_his)

  #igraph
  if(igraph.on==T){
  #igraph
  size.of.vertex<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length,na.action = stats::na.pass))
  #NA count as 1
  size.of.vertex1<-as.data.frame(stats::aggregate(barcode.history~seq.number,history,length))
  size.of.vertex<-merge(size.of.vertex1,size.of.vertex,all=T, by="seq.number")[,-3]
  size.of.vertex[!(size.of.vertex$seq.number %in% size.of.vertex1$seq.number),2]<-0
  #NA count as 0
  colnames(size.of.vertex)<-c("seq.number","size")
  rm(size.of.vertex1)

  #isotype distribution list
  history$isotype[is.na(history$barcode.history)]<-9999
  history$trans_state_his[is.na(history$barcode.history)]<-9999

  isotype.distribution<-reshape2::dcast(history, seq.number~isotype, value.var ="isotype",length)[,-1]
  if(length(which(colnames(isotype.distribution)==9999))!=0){
    isotype.distribution<-isotype.distribution[,1:ncol(isotype.distribution)-1]
    }
  isotype.name<-colnames(isotype.distribution)
  isotype.distribution <-as.list(as.data.frame(t(isotype.distribution)))
  #transcriptome state distribution list
  trans_state_distribution<-reshape2::dcast(history, seq.number~trans_state_his, value.var ="trans_state_his",length)[,-1]
  if(length(which(colnames(trans_state_distribution)==9999))!=0){
    trans_state_distribution<-trans_state_distribution[,1:ncol(trans_state_distribution)-1]
  }
  trans_state_name<-colnames(trans_state_distribution)
  trans_state_distribution <-as.list(as.data.frame(t(trans_state_distribution)))

  #write igraph list and assign attr
  for (i in 1:length(igraph.index)){

    if (length(igraph.index[[i]])>1){
      igraph.index[[i]]<-c(0,igraph.index[[i]])
      igraph.index.attr[[i]]<-data.frame(igraph.index[[i]])
      colnames(igraph.index.attr[[i]])<-"seq.number"

      igraph.index.jr[[i]]<-data.frame(igraph.index[[i]])
      colnames(igraph.index.jr[[i]])<-"seq.number"

      #delete duplicated for vertex attribute
      igraph.index.attr[[i]]<-subset(igraph.index.attr[[i]],!duplicated(igraph.index.attr[[i]]$seq.number))
      #size of vertex
      igraph.index.attr[[i]]<-merge(size.of.vertex,igraph.index.attr[[i]],by="seq.number",all.y=T)
      igraph.index.attr[[i]]$size[1]<-0
      igraph.index.attr[[i]]$adj_size<-(igraph.index.attr[[i]]$size*60/(range(igraph.index.attr[[i]]$size)[2]-range(igraph.index.attr[[i]]$size)[1]))^0.5+15 #visual size of vertex a number between 20 and 60
      #size of empty
      igraph.index.attr[[i]]$adj_size<-replace(igraph.index.attr[[i]]$adj_size, igraph.index.attr[[i]]$size==0, 10)

      #simplified index
      igraph.index.jr[[i]]$No<-match(x=igraph.index.jr[[i]]$seq.number,table=igraph.index.attr[[i]]$seq.number)

      #find empty
      size0<-which(igraph.index.attr[[i]]$size==0)

      No<- igraph.index.jr[[i]]$No

      #match isotype.distribution to each vertex
      pie.values<-list()
      for(j in 2:length(igraph.index.attr[[i]]$seq.number)){
        pie.values[[j]]<-isotype.distribution[[igraph.index.attr[[i]]$seq.number[j]]]
      }
      pie.values.list[[i]]<-pie.values

      #match trans state distribution to each vertex
      pie_values_trans<-list()
      for(j in 2:length(igraph.index.attr[[i]]$seq.number)){
        pie_values_trans[[j]]<-trans_state_distribution[[igraph.index.attr[[i]]$seq.number[j]]]
      }
      pie_values_trans_list[[i]]<-pie_values_trans

      #write graph list
      g<-igraph::graph(No,directed = T)
      g$layout<-igraph::layout_as_tree
      igraph::V(g)$label<-igraph.index.attr[[i]]$size
      igraph::V(g)$label[1]<-"Germline"
      igraph::V(g)$size<-igraph.index.attr[[i]]$adj_size
      igraph::V(g)$label.dist<-3

      igraph::V(g)$shape<-"pie"

      if(length(size0)>1){
        for(k in 1:length(size0)){
          igraph::V(g)$shape[[size0[k]]]<-"circle"
        }
      }
      if(length(size0)==1 ){
        igraph::V(g)$shape[[size0]]<-"circle"
      }

      igraph::V(g)$pie.color<-list(colors)
      igraph::V(g)$color<-"gray"
      igraph::V(g)$color[1]<-"black"
      p<-g
      ###!!!
      igraph::V(g)$pie<-pie.values.list[[i]]
      igraph::V(p)$pie<-pie_values_trans_list[[i]]



      igraph_list[[i]]<-g
      igraph_list_trans[[i]]<-p
    }

    if (length(igraph.index[[i]])==1){
      igraph.index.attr[[i]]<-size.of.vertex$size[which(size.of.vertex$seq.number==igraph.index[[i]])]#size
      igraph.index.jr[[i]]<-data.frame(igraph.index[[i]])
      pie.values.list[[i]]<-list(isotype.distribution[igraph.index[[i]]][[1]])
      pie_values_trans_list[[i]]<-list(trans_state_distribution[igraph.index[[i]]][[1]])

      g<-igraph::graph(c(1,1), edges = NULL)
      g$layout<-igraph::layout_as_tree
      igraph::V(g)$size<-igraph.index.attr[[i]]/30+20
      igraph::V(g)$label<-igraph.index.attr[[i]]
      igraph::V(g)$label.dist<-5
      igraph::V(g)$shape<-"pie"
      igraph::V(g)$pie.color<-list(colors)
      p<-g
      igraph::V(g)$pie<-pie.values.list[[i]]
      igraph::V(p)$pie<-pie_values_trans_list[[i]]

      igraph_list[[i]]<-g
      igraph_list_trans[[i]]<-p
    }
  }
  }
  if(igraph.on==F){
   igraph_list<-"none"
   igraph_list_trans<-"none"
 }

  #transcriptome
  if(transcriptome.on==T){
    ngene<-nrow(transcriptome.states)
    ncell<-length(barcode_uniq)
    transcriptome<-as.data.frame(matrix(0,ngene,ncell,dimnames = list(rownames(transcriptome.states),barcode_uniq)))
    for(i in 1:length(barcode_uniq)){
      transcriptome[,i]<-eval(parse(text=trans_dis[trans_state[i]]))
    }
    transcriptome<-as.matrix(transcriptome)
    transcriptome[transcriptome<0]<-0
  }
  if(transcriptome.on==F){
    transcriptome<-"none"
  }

  #Seq_Name

  Seq_Name<-"none"

  if(seq.name!=F){
    if(seq.name>clonotype_num){
    warning("seq.name > total number of clones")
    seq.name<-clonotype_num
    }
    c_gene_h<-c_gene[which(nrow(c_gene)%%2!=0)]
    Name<-paste0(clonotype_id,"_",c_gene_h,"_",barcode_uniq,"_cluster",trans_state)

    seq_germ<-c()
    for(i in 1:clonotype_num){
      seq_germ[i]<-paste(reference$reference_seq[2*i-1],reference$reference_seq[2*i],sep = "_")
    }

    clonotype_id_germ<-paste0("clonotype",1:clonotype_num)
    Name_germ<-rep("germline_clusterUnknown",clonotype_num)

    Seq<-c(seq_germ,seq_combi)
    Seq<-do::Replace(Seq,from = c(";IGK:",";IGL:"),to="_")
    Seq<-do::Replace(Seq,from = "IGH:",to="")

    Name<-c(Name_germ,Name)
    clonotype_id_temp<-c(clonotype_id_germ,clonotype_id)

    Seq_Name_temp<-data.frame(Seq,Name,clonotype_id_temp)
    frequency<-NULL
    uniq_clone<-dplyr::arrange(clonotypes,dplyr::desc(frequency))[1:seq.name,1]
    Seq_Name<-list()
  for(s in 1:length(uniq_clone)){
    Seq_Name[[s]]<- Seq_Name_temp[which(Seq_Name_temp$clonotype_id_temp==uniq_clone[s]),-3]
    }
  }

  #Version & Parameter
  version<-R.Version()
  parameters<-list(initial.size.of.repertoire,
  species,
  cell.type,
  cd4.proportion,
  duration.of.evolution,
  complete.duration,
  vdj.productive,
  vdj.model,
  vdj.insertion.mean,
  vdj.insertion.stdv,
  vdj.branch.prob,
  clonal.selection,
  cell.division.prob,
  sequence.selection.prob,
  special.v.gene,
  class.switch.prob,
  class.switch.selection.dependent,
  class.switch.independent,
  SHM.method,
  SHM.nuc.prob,
  SHM.isotype.dependent,
  SHM.phenotype.dependent,
  max.cell.number,
  max.clonotype.number,
  death.rate,
  igraph.on,
  transcriptome.on,
  transcriptome.switch.independent,
  transcriptome.switch.prob,
  transcriptome.switch.isotype.dependent,
  transcriptome.switch.SHM.dependent,
  transcriptome.switch.selection.dependent,
  transcriptome.states,
  transcriptome.noise,
  seq.name)

  names(parameters)<-c("initial.size.of.repertoire",
  "species",
  "cell.type",
  "cd4.proportion",
  "duration.of.evolution",
  "complete.duration",
  "vdj.productive",
  "vdj.model",
  "vdj.insertion.mean",
  "vdj.insertion.stdv",
  "vdj.branch.prob",
  "clonal.selection",
  "cell.division.prob",
  "sequence.selection.prob",
  "special.v.gene",
  "class.switch.prob",
  "class.switch.selection.dependent",
  "class.switch.independent",
  "SHM.method",
  "SHM.nuc.prob",
  "SHM.isotype.dependent",
  "SHM.phenotype.dependent",
  "max.cell.number",
  "max.clonotype.number",
  "death.rate",
  "igraph.on",
  "transcriptome.on",
  "transcriptome.switch.independent",
  "transcriptome.switch.prob",
  "transcriptome.switch.isotype.dependent",
  "transcriptome.switch.SHM.dependent",
  "transcriptome.switch.selection.dependent",
  "transcriptome.states",
  "transcriptome.noise",
  "seq.name")
#return----

  output_list<-list(all_contig_annotations,clonotypes,all_contig,consensus,reference,reference_real,transcriptome,igraph_list,igraph_list_trans,Seq_Name,igraph.index.attr,history,igraph.index,selected_seq,version,parameters)
  names(output_list)<-c("all_contig_annotations","clonotypes","all_contig","consensus","reference","reference_real","transcriptome","igraph_list_iso","igraph_list_trans","Seq_Name","igraph.index.attr","history","igraph.index","selected.seq","version","parameters")
  return(output_list)

}#function end
