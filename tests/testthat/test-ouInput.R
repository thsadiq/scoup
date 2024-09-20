# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

test_that("hbInput", {
    expect_is( ouInput(), "ou" )
    expect_error( ouInput("a") )
    expect_warning( ouInput(3) )
    expect_error( ouInput(c(eMean=100, eVar="10")) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #