# ><>< ================================================================ ><>< #
# ><><     scoup: Simulate Codon Sequences with Darwinian Selection     ><>< #
# ><><          Incorporated  as an Ornstein-Uhlenbeck Process          ><>< #
# ><>< ================================================================ ><>< #

test_that("aaGauss works properly", {
    expect_error( aaGauss(1, "") )
    expect_equal( unique( coeffs( aaGauss(0, 0) ) ), 0)
    expect_equal( var( coeffs( aaGauss(0.1, 1e-13) ) ), 0)
    expect_equal( all( coeffs( aaGauss(0.1, 0.2) ) >= 0 ), TRUE)
    expect_equal( unique( coeffs( aaGauss(0.1, 1e-13) ) ) > 0, TRUE)
    expect_equal( nsynVar( aaGauss(0.7, 1e-03) ) - 0.001, 0 )
})

# ><>< ================================================================ ><>< #
# ><><                          CODE ENDS HERE                          ><>< #
# ><>< ================================================================ ><>< #