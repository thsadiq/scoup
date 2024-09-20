# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

test_that("hbInput", {
    expect_error( hbInput("a") )
    expect_warning( hbInput(3) )
    expect_warning( hbInput(c(1, a=3,Ne=100)) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #