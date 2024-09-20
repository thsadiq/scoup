# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

# ><>< # Ornstein-Uhlenbeck Landscape Input
ouInput <- function(ouVector=0){
    stopifnot("`ouInput` arg isn't a numeric vector" = is(ouVector,"numeric"))
    warner(ouVector, c("eMean","eVar","Theta"))
    pull <- ouVector["eMean"]
    nPul <- ifelse(is.na(pull), 0, pull)
    vary <- ouVector["eVar"]
    nVar <- ifelse(is.na(vary), 0.01, vary)
    oscil <- ouVector["Theta"]
    nOsc <- ifelse(is.na(oscil), 0.01, oscil)
    scribe <- paste0("OU Parameters: Mu=",nPul,
                    "  Theta=",nOsc,"  Sigma.Sq=",nVar)
    names(nVar) <- NULL
    names(nOsc) <- NULL
    names(nPul) <- NULL
    ouEntry <- new("ou", var=nVar, theta=nOsc, mu=nPul, words=scribe)
    return(ouEntry)
}

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #