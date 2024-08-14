# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><      Amino Acid Selection Coefficient Generator (Gamma PDF)      ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 14 Aug, 2024                         ><>< #
# ><>< ================================================================ ><>< #

test_that("aaGamma works appropriately", {
    expect_error( aaGamma(1, "") )
    expect_equal( unique( aaGamma(0, 0) ), 0)
    expect_equal( var( aaGamma(0.1, 1e-13)[1:20] ), 0)
    expect_equal( all( aaGamma(0.1, 0.2)[1:20] >= 0 ), TRUE)
    expect_equal( unique( aaGamma(0.1, 1e-13)[1:20] ) > 0, TRUE)
    expect_equal( aaGamma(0.7, 1e-03)[22] - 0.001, c(nsynVar=0) )
    expect_equal( unique( names( aaGamma(0.1, 1e-13)[1:20] ) ), "")
    expect_equal( names( aaGamma(0.1, 1e-13)[21:22] ), c("synVar","nsynVar"))
})


# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #