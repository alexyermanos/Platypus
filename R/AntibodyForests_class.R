#' Class used for AntibodyForests functions
#'@description See AntibodyForests function for complete documentation
#' @export
#' @examples
#' \dontrun{
#' AntibodyForests(VDJ, sequence.type='VDJ.VJ.nt.trimmed',
#' include.germline=T, network.algorithm='tree',
#' resolve.ties=c('close.germline.distance', 'max.expansion'),
#' node.features='OVA_binder', expand.intermediates=T, network.level='intraclonal')
#'}

  methods::setClass('AntibodyForests',
    slots = c(
      tree = 'ANY', #in
      sample_id = 'ANY', #in
      clonotype_id = 'ANY', #in
      plot_ready = 'ANY', #f call
      heterogeneous = 'ANY', #f call
      reactivity = 'ANY', #f call
      dynamic = 'ANY', #f call
      metrics = 'ANY', #f call
      sequences = 'ANY', #in
      germline_sequence = 'ANY', #to add
      barcodes = 'ANY', #in
      node_features = 'ANY', #in
      edge_list = 'ANY', #no/f call
      gex_list = 'ANY', #f call
      paths = 'ANY', #f call
      node_transitions = 'ANY', #f call
      adjacency_matrix = 'ANY', #no/f call
      phylo = 'ANY', #f call
      feature_names = 'ANY',
      network_algorithm = 'ANY',
      inferred = 'ANY',
      permuted_transitions = 'ANY'
    )
  )

  methods::setMethod(f='show', signature='AntibodyForests',
   definition=function(object) {
    cat('AntibodyForests object', '\n')

    cat(length(which(object@node_features$node_type == 'sequence')), ' sequence nodes across ', length(object@sample_id), ' sample(s) and ', length(object@clonotype_id), ' clonotype(s)', '\n')

    if(any(object@node_features$node_type == 'intermediate')){
      cat(length(which(object@node_features$node_type == 'intermediate')), ' intermediate nodes', '\n')
    }

    if(any(object@node_features$node_type == 'inferred')){
      cat(length(which(object@node_features$node_type == 'inferred')), ' inferred sequence nodes', '\n')
    }

    if(any(object@node_features$node_type == 'bulk')){
      cat(length(which(object@node_features$node_type == 'bulk')), ' bulk sequence nodes', '\n')
    }

    cat(length(which(object@node_features$node_type == 'sequence')) + length(which(object@node_features$node_type == 'intermediate')) + length(which(object@node_features$node_type == 'inferred')) + length(which(object@node_features$node_type == 'bulk')),
        ' total network nodes', '\n'
    )


    cat('Sample id(s): ', paste0(object@sample_id, collapse = ', '), '\n')

    cat('Clonotype id(s): ', paste0(object@clonotype_id, collapse = ', '), '\n')

    networks <- c(object@network_algorithm)

    if(!is.null(object@plot_ready)){
      networks <- c(networks, 'plot_ready')
    }

    if(!is.null(object@phylo)){
      networks <- c(networks, 'phylo')
    }

    if(!is.null(object@heterogeneous)){
      networks <- c(networks, 'heterogeneous')

    }

    if(!is.null(object@reactivity)){
      networks <- c(networks, 'reactivity')

    }

    if(!is.null(object@dynamic)){
      networks <- c(networks, 'dynamic')
    }


    if(!is.null(object@metrics)){
      networks <- c(networks, 'metrics')

    }

    if(!is.null(object@node_transitions)){
      networks <- c(networks, 'node_transitions')

    }

    if(!is.null(object@paths)){
      networks <- c(networks, 'paths')
    }

    cat('Networks/analyses available: ', paste0(networks, collapse = ', '))

    cat('\n')

    cat('\n')

    cat('\n')

   }
  )
