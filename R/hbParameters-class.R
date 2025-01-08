# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # "hbParameters" Object Generics and Methods
setClass("hbParameters", representation(psize="numeric",
    vNvS="numeric", sampler="numeric", nsynVar="numeric", words="character"))

setMethod("show", "hbParameters", function(object){
    cat("\n", is(object)[[1]])
    cat("\nHalpern-Bruno Mutation-Selection Evolutionary Model Settings:")
    cat("\n\tPopulation size = ", object@psize)
    cat("\n\tNon-synonymous coefficients variance = ", object@nsynVar)
    cat("\n\tRatio of coefficients variances (vN/vS) = ", object@vNvS, "\n\n")
})
setMethod("vNvS", signature("hbParameters"), function(x) x@vNvS)
setMethod("effpop", signature("hbParameters"), function(x) x@psize)
setMethod("nsynVar", signature("hbParameters"), function(x) x@nsynVar)
setMethod("sampler", signature("hbParameters"), function(x) detectApp(x))

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #