#' Plot protein correlation network derived from coefficient matrix B
#'
#' @import tibble
#' @import magrittr
#' @import dplyr
#' @import Matrix
#' @import sna
#' @import graphics
#' @import network
#' @import ggnetwork
#' @import RColorBrewer
#' @export
#'
plot_network = function(tb_B,title_str) {

  # qunatiles over boostrap samples
  tb_B_summary = tb_B %>%
    group_by(protein_name) %>%
    summarize(
      coeff_median = median(coeff),
      min_median = quantile(coeff,probs = 0.05),
      max_median = quantile(coeff,probs = 0.95)
    )

  # extract nodes
  nodes_set1 = tb_B_summary %>% dplyr::filter(max_median < 0)
  nodes_set2 = tb_B_summary %>% dplyr::filter(min_median > 0)
  nodes_all = list(nodes_set1,nodes_set2) %>% bind_rows

  A = NULL
  if(nrow(nodes_all) > 0) {
    # make adjency matrices
    A1 = diag(nrow(nodes_set1))
    A1 = rbind(A1,rep(1,nrow(nodes_set1)))
    A1 = cbind(A1,rep(1,nrow(nodes_set1)+1))
    A2 = diag(nrow(nodes_set2))
    A2 = rbind(A2,rep(1,nrow(nodes_set2)))
    A2 = cbind(A2,rep(1,nrow(nodes_set2)+1))
    A = bdiag(A1,A2) %>% as.matrix
    pnames = c(nodes_set1$protein_name,"Hub1",
               nodes_set2$protein_name,"Hub2")
    dimnames(A) = list(pnames,pnames)
    super1 = which(colnames(A)=="Hub1")
    super2 = which(colnames(A)=="Hub2")
    A[super1,super2] = -1
    A[super2,super1] = -1
    if(sum(A[super1,]) == 0) A = A[-super1,-super1]
    super2 = which(colnames(A)=="Hub2")
    if(sum(A[super2,]) == 0) A = A[-super2,-super2]

    # make network
    net = network(A,
                  directed = FALSE,
                  ignore.eval = FALSE,
                  names.eval = "correlation")

    # ploting using ggnetwork (kamadakawai)
    net_df = ggnetwork(net,layout = "kamadakawai")
    net_df$correlation %<>% factor(levels = c(-1,1),
                                   labels = c("negative", "positive"))
    set.seed(0xdada2)
    remove_hubs = function(node) node[!node$vertex.names %in% c("Hub1","Hub2"), ]
    ggplot(net_df,aes(x = x,y = y,xend = xend,yend = yend)) +
      geom_edges(aes(color = correlation),size = 1,alpha = 0.6) +
      geom_nodes(size = 4,data = remove_hubs) +
      geom_nodelabel_repel(aes(label = vertex.names),data = remove_hubs) +
      #scale_color_brewer(palette = "Dark2") +
      scale_color_manual(values=c("#F8766D", "#00BFC4")) +
      theme_blank() +
      ggtitle(title_str)
  } else {
    ggplot() + ggtitle("no significant nodes found")
  }
}
