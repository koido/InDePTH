\name{cscore_LINCS}
\alias{cscore_LINCS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 Calculation CMap score from LINCS rank matrix.
}
\description{
 This function calculate CMap score by using LINCS rank matrix created from Landmark 978 genes. Before using this function, users must run init function only one time. Users set the output directory of init function as input directory of this function.
}
\usage{
cscore_LINCS(up.sig, dn.sig, input.dir, output.dir, write.name = NULL)
}
\arguments{
  \item{up.sig}{
  character vector of up sigunature of gene sets (Probe ID of Human Genome U133 Plus 2.0 Array is acceptable; if the arg 'gene.symbol' is set to TRUE, official gene symbol is also acceptable.).
}
  \item{dn.sig}{
  character vector of down sigunature of gene sets (Probe ID of Human Genome U133 Plus 2.0 Array is acceptable; if the arg 'gene.symbol' is set to TRUE, official gene symbol is also acceptable.).
}
  \item{input.dir}{
  a character string naming a directory containing the result of init function.
}
  \item{output.dir}{
  a character string naming a directory for output files from this function.
}
  \item{write.name}{
  a character string naming a file to write cscore results.
}
  \item{gene.symbol}{
  a boolean selecting whether "up.sig" and "dn.sig" are gene symbols or NOT.
  If TRUE, gene symbols are replaced to the probe id in the chip, because the rank matrix was made from probe level data (please see MATERIALS AND METHODS in our manuscript.).
  
}
}
\details{
}
\value{
}
\references{
  M. Koido, et al., 'Detection of hub genes for developing gene expression networks under anticancer drug treatment', Oncotarget, 9(49), 29097-29111, 2018.

  M. Ushijima, et al. 'Development of a gene expression database and related analysis programs for evaluation of anticancer compounds.' Cancer science 104(3), 360-368, 2013.

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
                     output.dir = "out",
                     write.name = NULL)
## End(Not run)
}
