#' Palette for SCISSOR plots
#'
#' @param set
#'
#' @import RColorBrewer wesanderson
#' @export
palette_SCISSOR = function(set="normal") {
  if (set=="normal") {
    palette(c(RColorBrewer::brewer.pal(8,"Set1"),
              RColorBrewer::brewer.pal(8,"Dark2"),
              wesanderson::wes_palette(n=4,"GrandBudapest1"),
              wesanderson::wes_palette(n=5,"Darjeeling1"),
              wesanderson::wes_palette(n=5,"Cavalcanti1")))
  } else if (set=="group") {
    palette(c(RColorBrewer::brewer.pal(5,"Set2"),
              RColorBrewer::brewer.pal(7,"Dark2"),
              RColorBrewer::brewer.pal(8,"Set3")[c(1,3:7)],
              wesanderson::wes_palette(n=4,"GrandBudapest1"),
              wesanderson::wes_palette(n=5,"Darjeeling1"),
              wesanderson::wes_palette(n=5,"Cavalcanti1")))
  }
}
