# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

test_that("wInput", {
    expect_error( wInput(3) )
    expect_warning( wInput(list(3)) )
    expect_error( wInput(list(vNvS=".1")) )
    expect_error( wInput(list(technique=3)) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #