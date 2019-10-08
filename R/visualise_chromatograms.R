#' Plot Extracted-ion chromatogram group.
#'
#' @importFrom tidyr gather
#' @importFrom ggplot2 ggplot ggtitle geom_vline geom_line theme theme_bw aes
#' @export
plotXICgroup <- function(XIC_group, peakAnnot = NULL, Title =NULL){
  df <- do.call("cbind", XIC_group)
  df <- df[,!duplicated(colnames(df))]
  df <- gather(df, key = "Transition", value = "Intensity", -time)
  g <- ggplot(df, aes(time, Intensity, col=Transition)) + geom_line(show.legend = FALSE) + theme_bw()
  if(!is.null(Title)) g <- g + ggtitle(paste0(Title)) + theme(plot.title = element_text(hjust = 0.5))
  if(!is.null(peakAnnot)){
    g <- g + geom_vline(xintercept=peakAnnot, lty="dotted", size = 0.4)
  }
  return(g)
}
