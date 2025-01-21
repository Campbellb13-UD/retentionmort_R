This package provides several avenues to predict and account for user-based mortality and tag loss during mark-recapture studies.
             When planning a study on a target species, the 
             `retentionmort_generation` function can be used to produce multiple
             synthetic mark-recapture datasets to anticipate the error 
             associated with a planned field study to guide method development
             to reduce error. Similarly, if field data was already collected, 
             the `retentionmort` function can be used to predict the error from
             already generated data to adjust for user-based mortality and tag
             loss. The `test_dataset_retentionmort` function will provide an
             example dataset of how data should be inputted into the function to
             run properly. Lastly, the `retentionmort_figure` function can be
             used on any dataset generated from either model function to produce
             an Rmarkdown prinout of preliminary analysis associated with the
             model, including summary statistics and figures.
