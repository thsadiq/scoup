# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

test_that("aaGamma works appropriately", {
    expect_error( aaGamma(1, "") )
    expect_equal( unique( coeffs( aaGamma(0, 0) ) ), 0)
    expect_equal( var( coeffs( aaGamma(0.1, 1e-13) ) ), 0)
    expect_equal( all( coeffs( aaGamma(0.1, 0.2) ) >= 0 ), TRUE)
    expect_equal( unique( coeffs( aaGamma(0.1, 1e-13) ) ) > 0, TRUE)
    expect_equal( nsynVar( aaGamma(0.7, 1e-03) ) - 0.001, 0 )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #