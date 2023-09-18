# Analysis in few words

To familiarize with analysis recommended to read vignette [here.](https://bcm-uga.github.io/TESS3_encho_sen/articles/main-vignette.html)

--------

Firstly, preparing the input data in form of transposed binary matrix and list of GPS coordinates from metadata for each isolate. Input files saved to RDS format (R specific)

### To instal tess3r pkg
```
install.packages("devtools")
devtools::install_github("bcm-uga/TESS3_encho_sen")
```

Additionally required packages: dplyr, readr, optparse, ggplot2, ggrepel

__Scripts:__

* `prepare_tess3r_files.R` (not parametrized) - prepare input RDS files. Required binary matrix, metadata with GPS coordinates

* `run_tess3r_files.R` - generate admxiture with tess3r for selected K ancestry

    * Single K

    ```
    Rscript ~/software/malaria-hub/tess3r/run_tess3r_files.R -b <binary_matrix.rds> -c <coords.rds> --prefix <prefix> --K 1 --rep 50 --threads 10
    ```

    * Multiple Ks

    ```
    cat Kruns.txt | xargs -I {} sh -c '~/software/malaria-hub/tess3r/run_tess3r_files.R -b <binary_matrix.rds> -c <coords.rds> --prefix <prefix> --K {} --rep 50 --threads 10'
    ```

* `plot_tess3r.R` (not parametrized) - generate admixture-like barplots and map
