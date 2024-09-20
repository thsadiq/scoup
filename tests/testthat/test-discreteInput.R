# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

test_that("discreteInput", {
    expect_error( discreteInput(NA) )
    expect_error( discreteInput(list(nodeIndex="Hi")) )
    expect_warning( discreteInput(list(qne="Hello", king=rep(1,3))) )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #