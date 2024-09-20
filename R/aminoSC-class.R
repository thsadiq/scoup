# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # "aminoSC" Object Generics and Methods
setClass("aminoSC",
    representation(coeffs="numeric", synVar="numeric", nsynVar="numeric"))

setMethod("coeffs", "aminoSC", function(x) nameDel(x@coeffs))
setMethod("synVar", "aminoSC", function(x) nameDel(x@synVar))
setMethod("nsynVar", "aminoSC", function(x) nameDel(x@nsynVar))
setMethod("show", "aminoSC", function(object) {
    cat("\n", is(object)[[1]])
    cat("\n\nAmino acid selection coefficients:\n")
    cat(paste0(aacids[seq(1,4)],
        "=",round(object@coeffs[seq(1,4)],3),", ")," ...\n")
    cat("\nSynonymous variance:\n",object@synVar,"\n")
    cat("\nNonsynonymous variance:\n",object@nsynVar,"\n\n")
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #