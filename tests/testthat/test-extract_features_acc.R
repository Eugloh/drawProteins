# trying to write some tests for my draw_proteins package
# my best schematic as of 20170630 only uses two functions
# extract_feat_acc
# and the phosphosite extraction function.

# the extract_feat_acc function works inside a loop which is NOT great.

# it works with multiple accession numbers and with one I think?


# write function to do unit tests with is extract_feat_acc
context("features_to_dataframe")

processed_protein_fixture <- function(sequence = paste(rep("A", 748),
                                                        collapse = "")) {
  list(
    accession = "TEST_PROCESSED",
    entryName = "TEST_PROCESSED_HUMAN",
    taxid = 9606,
    sequence = sequence,
    features = list(
      list(type = "SIGNAL", description = "Signal peptide",
           begin = 1, end = 19),
      list(type = "PROPEP", description = "Propeptide",
           begin = 20, end = 213),
      list(type = "CHAIN", description = "Mature chain",
           begin = 214, end = 748),
      list(type = "DOMAIN", description = "Peptidase M12B",
           begin = 220, end = 456),
      list(type = "TRANSMEM", description = "Helical",
           begin = 673, end = 693)
    )
  )
}

test_that("feature conversion preserves coordinates and sequence length", {
  converted <- feature_to_dataframe(list(processed_protein_fixture()))

  expect_type(converted$sequenceLength, "integer")
  expect_equal(converted$sequenceLength, rep(748L, 5))
  expect_equal(converted$type,
               c("SIGNAL", "PROPEP", "CHAIN", "DOMAIN", "TRANSMEM"))
  expect_equal(converted$begin, c(1, 20, 214, 220, 673))
  expect_equal(converted$end, c(19, 213, 748, 456, 693))
})

test_that("feature conversion records an unknown missing sequence length", {
  converted <- extract_feat_acc(processed_protein_fixture(sequence = NULL))

  expect_true(all(is.na(converted$sequenceLength)))
  expect_type(converted$sequenceLength, "integer")
})

test_that("missing CHAIN fallback spans the complete sequence inclusively", {
  fixture <- processed_protein_fixture(sequence = paste(rep("A", 100),
                                                        collapse = ""))
  fixture$features <- Filter(function(feature) feature$type != "CHAIN",
                             fixture$features)

  converted <- extract_feat_acc(fixture)
  chain <- converted[converted$type == "CHAIN", ]

  expect_equal(chain$begin, 1)
  expect_equal(chain$end, 100)
  expect_equal(chain$length, 100)
  expect_equal(chain$sequenceLength, 100L)
})

test_that("features_to_dataframe",{

  # load data from the package
  data("rel_json")
  # data for Rel A only
  rel_json %>%
    drawProteins::feature_to_dataframe() ->
    rel_data
  expect_is(rel_data, "data.frame")
  expect_equal(nrow(rel_data), 75)
  expect_equal(ncol(rel_data), 10)

  # load data from the package
  data("five_rel_list")
  # five_rel_list was created 20170818 using this code:
  # "Q04206 Q01201 Q04864 P19838 Q00653" %>%
  #   drawProteins::get_features() ->
  #   five_rel_list

  prot_data <- feature_to_dataframe(five_rel_list)

  expect_is(prot_data, "data.frame")

  # ncol at this point 10
  expect_equal(ncol(prot_data), 10)   # exact for sample data

  expect_equal(nrow(prot_data), 319)  # exact for sample data
  expect_match(prot_data[1,1], "CHAIN") # first element should be "CHAIN"
  expect_match(colnames(prot_data)[3], "begin") # 3rd column name is begin
  expect_match(colnames(prot_data)[4], "end") # 4th column name is end
  expect_match(colnames(prot_data)[5], "length") # 5th column name is length
  expect_match(colnames(prot_data)[6], "accession") # 6th col name - accession
  expect_match(colnames(prot_data)[7], "entryName") # 7th col name is entryName
  expect_match(colnames(prot_data)[8], "taxid") # 8th column name is taxid
  expect_match(colnames(prot_data)[9], "sequenceLength")
  expect_match(colnames(prot_data)[10], "order") # 10th column name is order

})


# write function to do unit tests with is extract_feat_acc
context("extract_feat_acc")

test_that("extract_feat_acc works to give ",{

  # load data from the package - should I move it into the test folder
  data("rel_A_features")  # this is object created from whole Uniprot GET
  # creates JSON object  in the environment.
  # it's the json data for Q04206 - for the transcription factor RelA

  # Nice simple test
  res <- extract_feat_acc(rel_A_features)

  expect_is(res, "data.frame") # generic

  # ncol at this point only 9 as order is added later.
  expect_equal(ncol(res), 9)   # exact for sample data

  expect_equal(nrow(res), 75)  # exact for sample data
  expect_match(res[1,1], "CHAIN") # first element should be "CHAIN"
  expect_match(colnames(res)[3], "begin") # 3rd column name is begin
  expect_match(colnames(res)[4], "end") # 4th column name is end
  expect_match(colnames(res)[6], "accession") # 6th column name is accession
  expect_match(colnames(res)[9], "sequenceLength")

  # need to add tests to check structure is going to work for plotting...

})



# write function to do unit tests with is extract_names
context("extract_names")

test_that("extract_names",{

  # load data from the package
  data("protein_json")
  prot_names <- extract_names(protein_json)

  expect_is(prot_names, "list")
  expect_equal(length(prot_names), 6)   # exact for sample data
  expect_equal(prot_names$accession, "Q04206")
  expect_equal(prot_names$name, "TF65_HUMAN")
  expect_equal(prot_names$protein.recommendedName.fullName,
                              "Transcription factor p65")
  expect_equal(prot_names$gene.name.primary, "RELA")
  expect_equal(prot_names$gene.name.synonym, "NFKB3")
  expect_equal(prot_names$organism.name.scientific, "Homo sapiens")
})
