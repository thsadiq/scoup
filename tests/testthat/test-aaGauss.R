# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><      Amino Acid Selection Coefficient Generator (Gauss PDF)      ><>< #
# ><><                       ~~~~~~~~~~~~~~~~~~~~                       ><>< #
# ><><                         V0: 14 Aug, 2024                         ><>< #
# ><>< ================================================================ ><>< #

test_that("aaGauss works properly", {
    expect_error( aaGauss(1, "") )
    expect_equal( unique( aaGauss(0, 0) ), 0)
    expect_equal( var( aaGauss(0.1, 1e-13)[1:20] ), 0)
    expect_equal( all( aaGauss(0.1, 0.2)[1:20] >= 0 ), TRUE)
    expect_equal( unique( aaGauss(0.1, 1e-13)[1:20] ) > 0, TRUE)
    expect_equal( aaGauss(0.7, 1e-03)[22] - 0.001, c(nsynVar=0) )
    expect_equal( unique( names( aaGauss(0.1, 1e-13)[1:20] ) ), "")
    expect_equal( names( aaGauss(0.1, 1e-13)[21:22] ), c("synVar","nsynVar"))
})


# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #