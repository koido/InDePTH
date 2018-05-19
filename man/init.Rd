\name{init_fun}
\alias{init_fun}
\title{
 init function of InDePTH.
}
\description{
      init_fun is ONLY requried for the first time of InDePTH analysis in your computer. It may use large memory and a lot of time (maybe more than a day). Following files must be set in the input.dir:

      zspc_n1328098x22268.gctx  /  level4 zspc data, provided from the lincscloud.org.

      q2norm_n1328098x22268.gctx  /  level3 q2norm data, provided from the lincscloud.org.

      inst.info  /  data description file, provided from the lincscloud.org.

      l000-chem-pert-20413 (1).xlsx  /  chemical information file obtained from http://support.lincscloud.org/hc/en-us/article_attachments/200662686/l000-chem-pert-20413__1_.xlsx

      HG-U133_Plus_2.na34.annot.csv  /  HG-U133_Plus_2 Annotations, CSV format, Release 34, obtained from Affymetrix.com support page.

}
\usage{
init(input.dir, out.dir)

}
\arguments{
  \item{input.dir}{
      a character string naming a directory containing input files.
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

  http://www.lincsproject.org/
}
\author{
  Masaru Koido (masaru.koido@jfcr.or.jp)
}
\examples{
## Not run:
init( input.dir = "input", out.dir = "init" )
## End(Not run)
}
