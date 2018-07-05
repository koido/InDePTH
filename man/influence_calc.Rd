\name{influence_calc}
\alias{influence_calc}
\title{
 InDePTH main function to calculate hub score.
}
\description{
 This function constructs DEGs network showing the degree of influence and from the network, calculate hub score. This function use the result object of cscore_calc function.
}
\usage{
 influence_calc(up.sig, dn.sig, up.ratio, dn.ratio, tot.thr, input.dir, output.dir, cscore = res, write.name = NULL)
}
\arguments{
  \item{up.sig}{
  character vector of up sigunature of gene sets (Probe ID of Human Genome U133 Plus 2.0 Array is acceptable).
}
  \item{dn.sig}{
  character vector of down sigunature of gene sets (Probe ID of Human Genome U133 Plus 2.0 Array is acceptable).
}
  \item{up.ratio}{
  numeric vector showing changes of gene expression levels of up sigunature of gene sets (probe id or official gene symbol).
}
  \item{dn.ratio}{
  numeric vector showing changes of gene expression levels of down sigunature of gene sets (probe id or official gene symbol).
}
  \item{tot.thr}{
  numeric of threshold to select genetic perturbation from the result of cscore_calc function.
}
  \item{input.dir}{
  a character string naming a directory containing the result of init function.
}
  \item{output.dir}{
  a character string naming a directory for output files from this function.
}
  \item{cscore}{
  a object from cscore_calc function.
}
  \item{write.name}{
  a character string naming a file to write network plot in igraph object (RDS format) and GML format. Suffix is NOT required.
}
}
\details{
}
\value{
}
\references{
  M. Koido, et al., 'Detection of hub genes for developing gene expression networks under anticancer drug treatment', Oncotarget, 9(49), 29097-29111, 2018.

  http://www.lincsproject.org/}
\author{
  Masaru Koido (masaru.koido@jfcr.or.jp)
}
\note{
}
\seealso{
}
\examples{
## Not run:
#-- only the first time of InDePTH
# init(input.dir = "input", out.dir = "init")
#
# Run (arbitrary threshold)
res1 <- cscore_LINCS(up.sig = c( "217761_at", "217398_x_at", "218744_s_at" ),
                     dn.sig = c( "221856_s_at", "218193_s_at", "218457_s_at" ),
                     input.dir = "init",
                     output.dir = "out",
                     write.name = NULL)
res2 <- influence_calc(up.sig = c( "217761_at", "217398_x_at", "218744_s_at" ),
                      dn.sig = c( "221856_s_at", "218193_s_at", "218457_s_at" ),
                      up.ratio = c( 2, 4, 6 ),
                      dn.ratio = c( 1/2, 1/3, 1/5 ),
                      tot.thr = 0.5,
                      input.dir = "init",
                      output.dir = "out",
                      cscore = res1,
                      write.name = NULL)

# Run (Data-driven threshold)
thr_res <- make_threshold( pert_id_vec = c( "TRCN0000010389", "TRCN0000010390", "TRCN0000010391"),
                           cscore = res1,
                           cell = "HT29",
                           input.dir = "init",
                           output.dir = "out" )
res2 <- influence_calc(up.sig = c( "217761_at", "217398_x_at", "218744_s_at" ),
                      dn.sig = c( "221856_s_at", "218193_s_at", "218457_s_at" ),
                      up.ratio = c( 2, 4, 6 ),
                      dn.ratio = c( 1/2, 1/3, 1/5 ),
                      tot.thr = thr_res$cutoff,
                      input.dir = "init",
                      output.dir = "out",
                      cscore = res1,
                      write.name = NULL)

## End(Not run)
}
