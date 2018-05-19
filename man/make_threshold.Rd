\name{make_threshold}
\alias{make_threshold}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Determine best cutoff point of CMap similarity score.
}
\description{
 This function calculate the best cutoff point of CMap score by using LINCS information file and user-determined true perturbation IDs. Before using this function, users must run init function only one time. Users set the output directory of init function as input directory of this function.
}
\usage{
make_threshold(pert_id_vec, cscore, cell, input.dir, output.dir)
}
\arguments{
  \item{pert_id_vec}{
  character vector of pert_id, described in LINCS inst.info file. LINCS data containing this pert_id vector are considered as true hit in ROC analysis.
}
  \item{cscore}{
  a object from cscore_calc function.
}
  \item{cell}{
  character vector of cell_line, described in LINCS inst.info file naming a directory containing the result of init function.
}
  \item{input.dir}{
  a character string naming a directory containing the result of init function.
}
  \item{output.dir}{
  a character string naming a directory for output files from this function.
}
}
\details{
}
\value{
}
\references{
  M. Koido, et al., 'Detection of hub genes for developing gene expression networks under anticancer drug treatment', under review.

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
# init( input.dir = "input", out.dir = "init" )
#
# Run
res1 <- cscore_LINCS(up.sig = c( "217761_at", "217398_x_at", "218744_s_at" ),
                     dn.sig = c( "221856_s_at", "218193_s_at", "218457_s_at" ),
                     input.dir = "init",
                     output.dir = "output",
                     write.name = NULL)
make_threshold( pert_id_vec = c( "TRCN0000010389", "TRCN0000010390", "TRCN0000010391"),
                cscore = res1,
                cell = "HT29",
                input.dir = "init",
                output.dir = "output" )
## End(Not run)
}
