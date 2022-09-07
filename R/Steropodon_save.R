Steropodon_save <- function(steropodon.object,
                            structure,
                            save.full,
                            save.dir
                           ){

   if(missing(steropodon.object)) stop('Please input your Steropodon object/ nested list of objects!')
   if(missing(structure)) structure <- 'structure'
   if(missing(save.full)) save.full <- T
   if(missing(save.dir)) save.dir <- 'steropodon_saved'


   if(!dir.exists(save.dir)){
     dir.create(save.dir)
   }

   save.dir <- paste0(save.dir, '/', structure)
   if(!dir.exists(save.dir)){
     dir.create(save.dir)
   }

   write_steropodon_pdbs(steropodon.object, structure = structure, dir = save.dir)

   if(save.full){
     steropodon_list <- unnest_steropodon(steropodon.object)
     file_names <- lapply(steropodon_list, function(x) paste0(x@structure_id, collapse = '_'))

     for(i in 1:length(steropodon_list)){
       struct <- select_structure(steropodon_list[[i]], structure = structure)
       file_path <- paste0(save.dir, '/', file_names[i], '.csv')
       write.csv(struct$atom, file_path)
     }
   }
}
