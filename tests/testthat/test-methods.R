context('methods test')

data("pbmc_small")

test_that("join_features",{


  pbmc_small %>% 
    join_features("CD3D") %>% 
    slice(1) %>%
    pull(.abundance_RNA) %>%
    expect_equal(6.35, tolerance=0.1)


})