#' Palette for SCISSOR plots
#'
#' @param set
#'
#' @import RColorBrewer wesanderson
#' @export
palette_SCISSOR = function(set="normal") {
  require(RColorBrewer); require(wesanderson)
  if (set=="normal") {
    palette(c(brewer.pal(8,"Set1"),brewer.pal(8,"Dark2"),
              wes_palette(n=4,"GrandBudapest1"),wes_palette(n=5,"Darjeeling1"),
              wes_palette(n=5,"Cavalcanti1")))
  } else if (set=="group") {
    palette(c(brewer.pal(5,"Set2"),brewer.pal(7,"Dark2"),brewer.pal(8,"Set3")[c(1,3:7)],
              wes_palette(n=4,"GrandBudapest1"),wes_palette(n=5,"Darjeeling1"),
              wes_palette(n=5,"Cavalcanti1")))
  }
}
