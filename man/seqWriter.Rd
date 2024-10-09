% File src/library/scoup/man/..Rd
\name{seqWriter}
\alias{seqWriter}

\title{Write Numeric Codon Alignment to a NEXUS File}
\description{
    Save numeric codon alignment matrix to a file in NEXUS format. It is
    particularly useful when data with site partitions is required.
}
\usage{
    seqWriter(alignmentMatrix, treeInfo=NA, addText="", fileTag=NULL)
}
\arguments{
    \item{alignmentMatrix}{A numerical matrix of codon sequence alignment
        that is similar to the \code{seqs} matrix from the output of
        \code{alignsim}. The rows of the matrix should each correspond to an
        extant taxa and the columns should be the alignments sites. The entries
        of the matrix should be integers between 1 and 61 and they will be
        decoded in terms of the ordered IUPAC sense codon triplets. That is,
        1=\code{AAA}, 2=\code{AAC}, 3=\code{AAG}, 4=\code{AAT}, 5=\code{ACA},
        \ldots, 57=\code{TGT}, 58=\code{TTA}, 59=\code{TTC}, 60=\code{TTG},
        \code{TTT}.}
    \item{treeInfo}{Phylogeny to be printed with the sequence. If
        unspecified (default = \code{NA}) a balanced phylogeny with branch
        length = 0.10 and number of extant taxa set as the number of rows of the
        input \code{alignmentMatrix} will be used.}
    \item{addText}{A string of comments to be printed with the alignment
        (default = "").}
    \item{fileTag}{Full path to where the output file should be printed. It
        should be a string (default = \code{NULL}). If not provided, the NEXUS
        file returned will be saved as \code{cranrSeqs.nex} in a temporary
        directory.}
}
\value{
    \describe{
        \item{\code{NULL}}{A NEXUS file with codon alignment printed therein
            will be saved in a temporary (or specified) directory.}
    }
}
\author{
    Hassan Sadiq
}
\seealso{
    Simulation function \code{\link{alignsim}}.
}
\examples{
    sqAlign <- alignsim(ouInput(), seqDetails(), hbInput(), NA)
    seqWriter(seqs(sqAlign))
}
