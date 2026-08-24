# unit tests for draw_chains
context("draw_canvas")
test_that("draw_canvas",{

  # load data from the package
  data("five_rel_data")
  # five_rel_data was created 20171101 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   five_rel_data

  # five_rel_data is a dataframe - 320 obs of 9 variables.
  p <- draw_canvas(five_rel_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have some labels
  expect_equal(p$labels$x, "Amino acid number")
  expect_equal(p$labels$y, "")

  # test default visual x limits
  expect_equal(p$coordinates$limits$x, c(-193.6, 1064.8))

  # test customizable visual x limits
  p_custom <- draw_canvas(five_rel_data, x_limits = c(-100, 1000))
  expect_equal(p_custom$coordinates$limits$x, c(-100, 1000))

})


# unit tests for draw_chains
context("draw_chains")

precursor_geometry_fixture <- function() {
  data.frame(
    type = c("CHAIN", "DOMAIN", "CHAIN", "CHAIN"),
    description = c("Mature chain", "Domain", "Mature chain", "Chain"),
    begin = c(214, 220, 214, 1),
    end = c(748, 456, 748, 100),
    length = c(534, 236, 534, 99),
    accession = c("TEST_PROCESSED", "TEST_PROCESSED",
                  "TEST_PROCESSED", "TEST_UNPROCESSED"),
    entryName = c("TEST_PROCESSED_HUMAN", "TEST_PROCESSED_HUMAN",
                  "TEST_PROCESSED_HUMAN", "TEST_UNPROCESSED_HUMAN"),
    taxid = 9606,
    sequenceLength = c(748L, 748L, 748L, 100L),
    order = c(1, 1, 2, 3),
    stringsAsFactors = FALSE
  )
}

test_that("precursor geometry creates one full-length segment per track", {
  result <- precursor_chain_data(precursor_geometry_fixture())

  expect_equal(nrow(result), 3)
  expect_equal(result$accession,
               c("TEST_PROCESSED", "TEST_PROCESSED", "TEST_UNPROCESSED"))
  expect_equal(result$order, c(1, 2, 3))
  expect_equal(result$begin, c(1, 1, 1))
  expect_equal(result$end, c(748, 748, 100))
})

test_that("precursor geometry rejects invalid or conflicting lengths", {
  data <- precursor_geometry_fixture()
  data$sequenceLength[2] <- 749L
  expect_error(precursor_chain_data(data),
               "TEST_PROCESSED.*order 1.*conflicting sequenceLength")

  for (bad_length in list(NA_real_, 0, -1, Inf)) {
    data <- precursor_geometry_fixture()
    data$sequenceLength[data$order == 3] <- bad_length
    expect_error(precursor_chain_data(data),
                 "TEST_UNPROCESSED.*order 3.*sequenceLength")
  }

  data <- precursor_geometry_fixture()
  data$sequenceLength <- NULL
  expect_error(precursor_chain_data(data),
               "sequenceLength.*regenerate.*feature_to_dataframe")
})

test_that("draw_chains",{

  # load data from the package
  data("five_rel_data")
  # five_rel_data was created 20171101 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   five_rel_data

  # five_rel_data is a dataframe - 320 obs of 9 variables.
  p <- draw_canvas(five_rel_data)
  p <- draw_chains(p, five_rel_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have two layers at this point
  expect_equal(length(p$layers), 2)
# https://stackoverflow.com/questions/13457562/how-to-determine-the-geom-type-of-each-layer-of-a-ggplot2-object/43982598#43982598
  # types of layers...
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomSegment")  # drawprotein chain here
  expect_equal(class(p$layers[[2]]$geom)[1], "GeomText")
  expect_equal(length(five_rel_data[five_rel_data$type == "DOMAIN",]),
              length(p$layers[[1]]$data))
})

test_that("draw_chains defaults to annotated chains and opts into precursor", {
  data <- precursor_geometry_fixture()
  original <- data

  annotated <- draw_chains(draw_canvas(data), data, label_chains = FALSE)
  expect_equal(annotated$layers[[1]]$data$begin, c(214, 214, 1))
  expect_equal(annotated$layers[[1]]$data$end, c(748, 748, 100))

  precursor <- draw_chains(draw_canvas(data), data,
                           extent = "precursor", label_chains = FALSE)
  expect_equal(precursor$layers[[1]]$data$order, c(1, 2, 3))
  expect_equal(precursor$layers[[1]]$data$begin, c(1, 1, 1))
  expect_equal(precursor$layers[[1]]$data$end, c(748, 748, 100))
  expect_equal(data, original)
  expect_equal(length(precursor$layers), 1)
})

test_that("draw_chains validates extent", {
  data <- precursor_geometry_fixture()
  expect_error(draw_chains(draw_canvas(data), data, extent = "complete"),
               "arg.*annotated.*precursor")
})



# unit tests with draw_domains
context("draw_domains")

test_that("draw_domains",{

  # load data from the package
  data("five_rel_data")
  # five_rel_data was created 20171101 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   five_rel_data

  # five_rel_data is a dataframe - 320 obs of 9 variables.
  p <- draw_canvas(five_rel_data)
  p <- draw_chains(p, five_rel_data)
  p <- draw_domains(p, five_rel_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have four layers at this point
  expect_equal(length(p$layers), 4)
  # two from draw_chains and two from draw_domains

  expect_equal(length(five_rel_data[five_rel_data$type == "DOMAIN",]),
              length(p$layers[[3]]$data))

  # types of layers, 3 and 4 added by draw_domains
  # expect_equal(class(p$layers[[3]]$geom)[1], "GeomRect")
  # expect_equal(class(p$layers[[4]]$geom)[1], "GeomLabel")
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomInteractiveRect") # ggiraph
  expect_equal(class(p$layers[[4]]$geom)[1], "GeomInteractiveLabel") # ggiraph
})

test_that("draw_domains renders precursor annotations by requested type", {
  data <- data.frame(
    type = c("SIGNAL", "PROPEP", "DOMAIN"),
    description = c("Signal peptide", "Propeptide", "Peptidase M12B"),
    begin = c(1, 20, 220),
    end = c(19, 213, 456),
    order = 1,
    stringsAsFactors = FALSE
  )
  canvas <- ggplot2::ggplot()

  signal <- draw_domains(canvas, data, type = "SIGNAL",
                         label_domains = FALSE)
  expect_equal(signal$layers[[1]]$data$type, "SIGNAL")
  expect_equal(signal$layers[[1]]$data$begin, 1)
  expect_equal(signal$layers[[1]]$data$end, 19)
  # expect_equal(rlang::as_label(signal$layers[[1]]$mapping$tooltip),
  #              "description") ##  importations '::' ou ':::' non déclarées depuis : ‘rlang’

  propeptide <- draw_domains(canvas, data, type = "PROPEP",
                             label_domains = FALSE)
  expect_equal(propeptide$layers[[1]]$data$type, "PROPEP")
  expect_equal(propeptide$layers[[1]]$data$begin, 20)
  expect_equal(propeptide$layers[[1]]$data$end, 213)

  absent <- draw_domains(canvas, data, type = "MOTIF",
                         label_domains = FALSE)
  expect_equal(nrow(absent$layers[[1]]$data), 0)
})




# unit tests with draw_phospho
context("draw_phospho")

test_that("draw_phospho",{

  # load data from the package
  data("five_rel_data")
  # five_rel_data was created 20171101 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   five_rel_data

  # five_rel_data is a dataframe - 320 obs of 9 variables.
  p <- draw_canvas(five_rel_data)
  p <- draw_chains(p, five_rel_data)
  p <- draw_phospho(p, five_rel_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have three layers at this point
  expect_equal(length(p$layers), 3)
  # two from draw_chains and one from draw_phospho
  # layers, 3 and 4 added by draw_domains
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomPoint")

  # should be 32 phosphorylation sites across data set...
  expect_equal(32, nrow(p$layers[[3]]$data))

  })

# unit tests with draw_motif
context("draw_motif")

test_that("draw_motif",{

  # load data from the package
  data("five_rel_data")
  # five_rel_data was created 20171101 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   five_rel_data

  # five_rel_data is a dataframe - 320 obs of 9 variables.
  p <- draw_canvas(five_rel_data)
  p <- draw_chains(p, five_rel_data)
  p <- draw_motif(p, five_rel_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have three layers at this point
  expect_equal(length(p$layers), 3)
  # two from draw_chains and one from draw_motif
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomRect")
  # https://stackoverflow.com/questions/13457562/how-to-determine-the-geom-type-of-each-layer-of-a-ggplot2-object/43982598#43982598
  # types of layers...
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomSegment")  # drawprotein chain here
  expect_equal(class(p$layers[[2]]$geom)[1], "GeomText")

  # p should have some labels
  expect_equal(p$labels$x, "Amino acid number")
  expect_equal(p$labels$y, "")
  expect_equal(length(five_rel_data[five_rel_data$type == "DOMAIN",]),
              length(p$layers[[1]]$data))
  # p$layers[[3]]$data contains the data that was extracted
  # dimensions are 6 9
  expect_equal(nrow(five_rel_data[five_rel_data$type == "MOTIF",]),
              nrow(p$layers[[3]]$data))
  expect_equal(p$layers[[3]]$data$type[1], "MOTIF" )
})


# unit tests with draw_regions
context("draw_regions")

test_that("draw_regions",{

  # load data from the package
  data("five_rel_data")
  # five_rel_data was created 20171101 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   five_rel_data

  # five_rel_data is a dataframe - 320 obs of 9 variables.
  p <- draw_canvas(five_rel_data)
  p <- draw_chains(p, five_rel_data)
  p <- draw_regions(p, five_rel_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have three layers at this point
  expect_equal(length(p$layers), 3)
  # two from draw_chains and one from draw_regions
  # https://stackoverflow.com/questions/13457562/how-to-determine-the-geom-type-of-each-layer-of-a-ggplot2-object/43982598#43982598
  # types of layers...
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomSegment")  # drawprotein chain here
  expect_equal(class(p$layers[[2]]$geom)[1], "GeomText")
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomRect")

  # p should have some labels
  expect_equal(p$labels$x, "Amino acid number")
  expect_equal(p$labels$y, "")
  expect_equal(length(five_rel_data[five_rel_data$type == "DOMAIN",]),
              length(p$layers[[1]]$data))
  # p$layers[[3]]$data contains the data that was extracted
  # dimensions are 6 9
  expect_equal(nrow(five_rel_data[five_rel_data$type == "REGION",]),
              nrow(p$layers[[3]]$data))
  expect_equal(p$layers[[3]]$data$type[1], "REGION" )
})



# unit tests with draw_repeat
context("draw_repeat")

test_that("draw_repeat",{

  # load data from the package
  data("five_rel_data")
  # five_rel_data was created 20171101 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   five_rel_data

  # five_rel_data is a dataframe - 320 obs of 9 variables.
  p <- draw_canvas(five_rel_data)
  p <- draw_chains(p, five_rel_data)
  p <- draw_repeat(p, five_rel_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have four layers at this point
  expect_equal(length(p$layers), 4)
  # two from draw_chains and two from draw_repeat
  # first draw_repeat layer is rectanges
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomRect")
  # second draw_repeat layer is text
  expect_equal(class(p$layers[[4]]$geom)[1], "GeomText")

  # p should have some labels
  expect_equal(p$labels$x, "Amino acid number")
  expect_equal(p$labels$y, "")
  expect_equal(length(five_rel_data[five_rel_data$type == "REPEAT",]),
              length(p$layers[[1]]$data))
  # p$layers[[3]]$data contains the data that was extracted
  # dimensions are 6 9
  expect_equal(nrow(five_rel_data[five_rel_data$type == "REPEAT",]),
              nrow(p$layers[[3]]$data))
  expect_equal(p$layers[[3]]$data$type[1], "REPEAT" )
})


# unit tests for draw_recept_dom
context("draw_recept_dom")

test_that("draw_recept_dom",{

  # load data from the package
  data("tnfs_data")
  # tnfs_data was created 20190103 using this code:
  # "P19438 P25942" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   tnfs_data

  # tnfs_data is a dataframe - 127 obs of 9 variables.
  p <- draw_canvas(tnfs_data)
  p <- draw_chains(p, tnfs_data)
  p <- draw_recept_dom(p, tnfs_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have four layers at this point
  # because TOPO_DOM and TRANSMEM added separately
  expect_equal(length(p$layers), 4)
  # https://stackoverflow.com/questions/13457562/how-to-determine-the-geom-type-of-each-layer-of-a-ggplot2-object/43982598#43982598
  # types of layers...
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomSegment")  # drawprotein chain here
  expect_equal(class(p$layers[[2]]$geom)[1], "GeomText")
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomRect")
  expect_equal(class(p$layers[[4]]$geom)[1], "GeomRect")

  # in this case nrow of layer 3 should be 4
  # and nrow of layer 4 should be 2
  expect_equal(nrow(p$layers[[3]]$data), 4) # 2 receptors - EC and Cyto domain
  expect_equal(nrow(p$layers[[4]]$data), 2) # 2 receptors - 2 TM domain


  ## tests for label domains (default is FALSE)
  # tnfs_data is a dataframe - 127 obs of 9 variables.
  p <- draw_canvas(tnfs_data)
  p <- draw_chains(p, tnfs_data)
  p <- draw_recept_dom(p, tnfs_data, label_domains = TRUE)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have six layers at this point
  # because TOPO_DOM and TRANSMEM added separately
  # and labelling domains gives two more layers which are GeomLable
  expect_equal(length(p$layers), 6)
  # https://stackoverflow.com/questions/13457562/how-to-determine-the-geom-type-of-each-layer-of-a-ggplot2-object/43982598#43982598
  # types of layers...
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomSegment")  # drawprotein chain here
  expect_equal(class(p$layers[[2]]$geom)[1], "GeomText")
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomRect")
  expect_equal(class(p$layers[[4]]$geom)[1], "GeomRect")
  expect_equal(class(p$layers[[5]]$geom)[1], "GeomLabel")
  expect_equal(class(p$layers[[6]]$geom)[1], "GeomLabel")

  # in this case nrow of layer 5 should be 2
  # and nrow of layer 6 should be 2
  expect_equal(nrow(p$layers[[5]]$data), 4) # 2 receptors - EC and Cyto domain
  expect_equal(nrow(p$layers[[6]]$data), 2) # 2 receptors - 2 TM domain

  # testing data and mapping of layers 5 and 6
  expect_equal(length(p$layers[[5]]$data), 9)
  expect_equal(nrow(p$layers[[5]]$data), 4)
  expect_equal(p$layers[[5]]$data$type[1], "TOPO_DOM")
  expect_equal(p$layers[[5]]$data$description[1], "Extracellular")
  expect_equal(p$layers[[5]]$data$description[2], "Cytoplasmic")
  expect_equal(p$layers[[6]]$mapping$label, "TM")
  expect_equal(p$layers[[6]]$data$type[1], "TRANSMEM")
  expect_equal(p$layers[[6]]$data$description[1], "Helical")

})


# useful advice here:
#https://stackoverflow.com/questions/31038709/how-to-write-a-test-for-a-ggplot-plot


# unit tests for draw_folding
context("draw_folding")

test_that("draw_folding",{

  # load data from the package
  data("tnfs_data")
  # tnfs_data was created 20190103 using this code:
  # "P19438 P25942" %>%
  #   drawProteins::get_features() %>%
  #   drawProteins::feature_to_dataframe() ->
  #   tnfs_data

  # tnfs_data is a dataframe - 127 obs of 9 variables.
  p <- draw_canvas(tnfs_data)
  p <- draw_chains(p, tnfs_data)
  p <- draw_folding(p, tnfs_data)

  # p is a ggplot object
  expect_is(p,"ggplot")

  # p should have five layers at this point
  # because HELIX, STRAND and TURN added separately
  expect_equal(length(p$layers), 5)
  # https://stackoverflow.com/questions/13457562/how-to-determine-the-geom-type-of-each-layer-of-a-ggplot2-object/43982598#43982598
  # types of layers...
  expect_equal(class(p$layers[[1]]$geom)[1], "GeomSegment")  # drawprotein chain here
  expect_equal(class(p$layers[[2]]$geom)[1], "GeomText")
  expect_equal(class(p$layers[[3]]$geom)[1], "GeomRect")
  expect_equal(class(p$layers[[4]]$geom)[1], "GeomRect")
  expect_equal(class(p$layers[[5]]$geom)[1], "GeomRect")

  # in this case nrow of layer 3 should be 4
  # and nrow of layer 4 should be 2
  expect_equal(nrow(p$layers[[3]]$data), 34) # 34 STRAND regions
  expect_equal(nrow(p$layers[[4]]$data), 10) # 10 HELIX regions
  expect_equal(nrow(p$layers[[5]]$data), 1) # 1 TURN

})
