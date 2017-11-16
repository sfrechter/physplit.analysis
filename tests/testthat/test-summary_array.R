context('summary array')

test_that("create_summary_array works",{
  cells=c("nm20130208c1","nm20130211c3")
  expect_is(res <- create_raw_summary_array(cells = cells), 'array')
  dn <-
    list(
      c("nm20130208c1", "nm20130211c3"),
      c(
        "OilBl",
        "E2Hex",
        "GerAc",
        "Prpyl",
        "IPenA",
        "Et3HB",
        "Nonnl",
        "CiVAc",
        "MetSl",
        "HexAc",
        "PeEtA",
        "AceAc",
        "EtHex",
        "2PnAc",
        "5OdMx",
        "BeZal",
        "bCitr",
        "1HxOl",
        "Frnsl",
        "WatBl",
        "Cdvrn",
        "Sprmn",
        "Acoin",
        "MtAct",
        "AcAcd",
        "PrpnA",
        "BtrAc",
        "Amnia",
        "Pyrdn",
        "PAcHd",
        "HCL36",
        "PAcAc",
        "Vingr",
        "Geosn",
        "VinGe",
        "PEtAm",
        "ClrBL",
        "ClrB2",
        "FlyFM",
        "EtAmn",
        "MtAmn",
        "Ptscn",
        "Lnlol",
        "23BTD"
      ),
      c("baseline",
        "max1", "max2", "max3", "max4", "max5", "max6")
    )
  expect_equal(dimnames(res), dn)
})
