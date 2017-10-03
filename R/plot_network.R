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
    A1 = matrix(data = 1,
                nrow = nrow(nodes_set1),
                ncol = nrow(nodes_set1))
    A2 = matrix(data = 1,
                nrow = nrow(nodes_set2),
                ncol = nrow(nodes_set2))
    A = bdiag(A1,A2) %>% as.matrix
    A[A==0] = -1
    dimnames(A) = list(nodes_all$protein_name,nodes_all$protein_name)
    
    # make network
    net = network(A,
                  directed = FALSE,
                  ignore.eval = FALSE,
                  names.eval = "correlation")
    
    # ploting using ggnetwork
    net_df = ggnetwork(net,layout = "circle")
    net_df$correlation %<>% factor(levels = c(-1,1),
                                   labels = c("negative", "positive"))
    set.seed(0xdada2)
    ggplot(net_df,aes(x = x,y = y,xend = xend,yend = yend)) +
      geom_edges(aes(color = correlation),size = 1,alpha = 0.6) +
      geom_nodes(size = 4) +
      geom_nodelabel_repel(aes(label = vertex.names)) +
      #scale_color_brewer(palette = "Dark2") +
      scale_color_manual(values=c("#F8766D", "#00BFC4")) +
      theme_blank() +
      ggtitle(title_str) +
      coord_fixed(ratio = 0.8)
  } else {
    ggplot() + ggtitle("no significant nodes found")
  }
}
