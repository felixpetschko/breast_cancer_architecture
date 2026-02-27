``` r
library(spacedeconv)
library(ggplot2)
library(cowplot)

dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

slides <- c("andersson", "1142243F", "1160920F", "4465", "44971", "4290", "4535")
title_size <- 18
font_size <- 14
legend_size <- 14
```

## Methods

Normalization and deconvolution follow the spacedeconv paper scripts:

- `preprocess(min_umi = 87)`
- `spacedeconv::normalize(...)`
- `deconvolute(method = "estimate", assay_sp = "cpm")`
- `deconvolute(method = "quantiseq", assay_sp = "cpm", tumor = TRUE)`
- `deconvolute(method = "epic", assay_sp = "cpm", tumor = TRUE)`

## Per-slide Spatial Results

``` r
for (slide in slides) {
  cat("\n### Slide ", slide, "\n\n", sep = "")

  est <- readRDS(file.path("results/objects", paste0("deconv_estimate_", slide, ".rds")))
  qua <- readRDS(file.path("results/objects", paste0("deconv_quantiseq_", slide, ".rds")))
  epi <- readRDS(file.path("results/objects", paste0("deconv_epic_", slide, ".rds")))

  p_tumor <- plot_spatial(
    est,
    result = "estimate_tumor.purity",
    smooth = TRUE,
    density = FALSE,
    title = paste0(slide, " - Tumor Purity (ESTIMATE)"),
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
  )

  p_cd8 <- plot_spatial(
    qua,
    result = "quantiseq_T.cell.CD8.",
    smooth = TRUE,
    density = FALSE,
    shift_positive = FALSE,
    title = paste0(slide, " - CD8 T cells (quanTIseq)"),
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
  )

  p_caf <- plot_spatial(
    epi,
    result = "epic_Cancer.associated.fibroblast",
    smooth = TRUE,
    density = FALSE,
    title = paste0(slide, " - CAFs (EPIC)"),
    title_size = title_size,
    font_size = font_size,
    legend_size = legend_size
  )

  panel <- cowplot::plot_grid(p_tumor, p_cd8, p_caf, ncol = 3, align = "hv")
  print(panel)

  ggsave(
    filename = file.path("results/plots", paste0("slide_", slide, "_targets.png")),
    plot = panel,
    width = 16,
    height = 5,
    units = "in",
    dpi = 300,
    bg = "white"
  )
}
```

    ## 
    ## ### Slide andersson

![](/gpfs/gpfs1/scratch/c9881013/felix_breast_cancer/breast_cancer/report/breast_spacedeconv_analysis_files/figure-gfm/per-slide-plots-1.png)<!-- -->

    ## 
    ## ### Slide 1142243F

![](/gpfs/gpfs1/scratch/c9881013/felix_breast_cancer/breast_cancer/report/breast_spacedeconv_analysis_files/figure-gfm/per-slide-plots-2.png)<!-- -->

    ## 
    ## ### Slide 1160920F

![](/gpfs/gpfs1/scratch/c9881013/felix_breast_cancer/breast_cancer/report/breast_spacedeconv_analysis_files/figure-gfm/per-slide-plots-3.png)<!-- -->

    ## 
    ## ### Slide 4465

![](/gpfs/gpfs1/scratch/c9881013/felix_breast_cancer/breast_cancer/report/breast_spacedeconv_analysis_files/figure-gfm/per-slide-plots-4.png)<!-- -->

    ## 
    ## ### Slide 44971

![](/gpfs/gpfs1/scratch/c9881013/felix_breast_cancer/breast_cancer/report/breast_spacedeconv_analysis_files/figure-gfm/per-slide-plots-5.png)<!-- -->

    ## 
    ## ### Slide 4290

![](/gpfs/gpfs1/scratch/c9881013/felix_breast_cancer/breast_cancer/report/breast_spacedeconv_analysis_files/figure-gfm/per-slide-plots-6.png)<!-- -->

    ## 
    ## ### Slide 4535

![](/gpfs/gpfs1/scratch/c9881013/felix_breast_cancer/breast_cancer/report/breast_spacedeconv_analysis_files/figure-gfm/per-slide-plots-7.png)<!-- -->

## Combined Table Preview

``` r
tab <- read.csv("results/tables/targets_all_slides.csv")
knitr::kable(head(tab, 20))
```

| spot_id | slide_id | tumor_purity_estimate | cd8_quantiseq | caf_epic | pxl_col_in_fullres | pxl_row_in_fullres |
|:---|:---|---:|---:|---:|---:|---:|
| AAACAAGTATCTCCCA-1 | andersson | 0.5339460 | 0.0871836 | 0.2495293 | 17428 | 15937 |
| AAACACCAATAACTGC-1 | andersson | 0.8913474 | 0.0000000 | 0.0419921 | 6092 | 18054 |
| AAACAGAGCGACTCCT-1 | andersson | 0.7072153 | 0.0000000 | 0.1191177 | 16351 | 7383 |
| AAACAGGGTCTATATT-1 | andersson | 0.5154692 | 0.1148961 | 0.0424538 | 5278 | 15202 |
| AAACAGTGTTCCTGGG-1 | andersson | 0.8476866 | 0.0000000 | 0.1106762 | 9363 | 21386 |
| AAACATTTCCCGGATT-1 | andersson | 0.8841088 | 0.0324866 | 0.0301322 | 16740 | 18549 |
| AAACCCGAACGAAATC-1 | andersson | 0.5437192 | 0.0000000 | 0.2154445 | 19205 | 14752 |
| AAACCGGGTAGGTACC-1 | andersson | 0.9024869 | 0.0000000 | 0.0238491 | 7328 | 14018 |
| AAACCTAAGCAGCCGG-1 | andersson | 0.8391736 | 0.0000000 | 0.0974704 | 14827 | 19495 |
| AAACCTCATGAAGTTG-1 | andersson | 0.7184005 | 0.0615161 | 0.1841585 | 6102 | 12828 |
| AAACGAAGAACATACC-1 | andersson | 0.8814044 | 0.0688613 | 0.0241004 | 12259 | 5475 |
| AAACGAGACGGTTGAT-1 | andersson | 0.7612716 | 0.0275129 | 0.0913638 | 14294 | 12368 |
| AAACGCCCGAGATCGG-1 | andersson | 0.7687589 | 0.0000000 | 0.0516603 | 18267 | 5011 |
| AAACGGGCGTACGGGT-1 | andersson | 0.8917116 | 0.0182992 | 0.0301744 | 15919 | 19497 |
| AAACGGTTGCGAACTG-1 | andersson | 0.7230822 | 0.0000000 | 0.1443731 | 11550 | 19964 |
| AAACGTGTTCGCCCTA-1 | andersson | 0.6892423 | 0.0000000 | 0.0681070 | 19628 | 7389 |
| AAACTAACGTGGCGAC-1 | andersson | 0.7234872 | 0.0000000 | 0.0995721 | 18538 | 5962 |
| AAACTCGGTTCGCAAT-1 | andersson | 0.9032156 | 0.0000000 | 0.0476065 | 13052 | 19730 |
| AAACTCGTGATATAAG-1 | andersson | 0.7274750 | 0.0000000 | 0.2223779 | 18941 | 9526 |
| AAACTGCTGGCTCCAA-1 | andersson | 0.7585878 | 0.0000000 | 0.0612951 | 12651 | 14740 |

## Session Info

``` r
sessionInfo()
```

    ## R version 4.3.3 (2024-02-29)
    ## Platform: x86_64-conda-linux-gnu (64-bit)
    ## Running under: Rocky Linux 8.10 (Green Obsidian)
    ## 
    ## Matrix products: default
    ## BLAS/LAPACK: /gpfs/gpfs1/scratch/c9881013/.conda_envs/spacedeconv-env/lib/libmkl_rt.so;  LAPACK version 3.8.0
    ## 
    ## locale:
    ##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
    ##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
    ##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
    ##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
    ##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    ## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    ## 
    ## time zone: Europe/Vienna
    ## tzcode source: system (glibc)
    ## 
    ## attached base packages:
    ## [1] stats4    stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] cowplot_1.2.0               ggplot2_3.5.2              
    ##  [3] rmarkdown_2.29              EPIC_1.1.7                 
    ##  [5] SpatialExperiment_1.12.0    SingleCellExperiment_1.24.0
    ##  [7] SummarizedExperiment_1.32.0 Biobase_2.62.0             
    ##  [9] GenomicRanges_1.54.1        GenomeInfoDb_1.38.1        
    ## [11] IRanges_2.36.0              S4Vectors_0.40.2           
    ## [13] BiocGenerics_0.48.1         MatrixGenerics_1.14.0      
    ## [15] matrixStats_1.5.0           spacedeconv_1.0.0          
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] splines_4.3.3           later_1.4.4             bitops_1.0-9           
    ##   [4] filelock_1.0.3          tibble_3.3.0            cellranger_1.1.0       
    ##   [7] preprocessCore_1.64.0   XML_3.99-0.17           lifecycle_1.0.4        
    ##  [10] sf_1.0-20               rstatix_0.7.2           edgeR_4.0.16           
    ##  [13] lattice_0.22-7          MASS_7.3-60.0.1         backports_1.5.0        
    ##  [16] magrittr_2.0.3          limma_3.58.1            yaml_2.3.10            
    ##  [19] remotes_2.5.0           reticulate_1.43.0       DBI_1.2.3              
    ##  [22] RColorBrewer_1.1-3      lubridate_1.9.4         abind_1.4-5            
    ##  [25] zlibbioc_1.48.0         rvest_1.0.5             quadprog_1.5-8         
    ##  [28] purrr_1.1.0             RCurl_1.98-1.17         pracma_2.4.4           
    ##  [31] rappdirs_0.3.3          sva_3.50.0              multimode_1.5          
    ##  [34] circlize_0.4.16         GenomeInfoDbData_1.2.11 data.tree_1.2.0        
    ##  [37] terra_1.8-42            genefilter_1.84.0       units_0.8-7            
    ##  [40] annotate_1.80.0         codetools_0.2-20        DelayedArray_0.28.0    
    ##  [43] xml2_1.4.0              tidyselect_1.2.1        shape_1.4.6.1          
    ##  [46] farver_2.1.2            BiocFileCache_2.10.1    jsonlite_2.0.0         
    ##  [49] e1071_1.7-16            ks_1.15.1               Formula_1.2-5          
    ##  [52] ggridges_0.5.7          survival_3.8-3          systemfonts_1.2.3      
    ##  [55] quantiseqr_1.10.0       tools_4.3.3             progress_1.2.3         
    ##  [58] ragg_1.5.0              Rcpp_1.1.0              glue_1.8.0             
    ##  [61] mnormt_2.1.1            SparseArray_1.2.2       xfun_0.53              
    ##  [64] mgcv_1.9-3              decoupleR_2.8.0         dplyr_1.1.4            
    ##  [67] withr_3.0.2             ComICS_1.0.4            fastmap_1.2.0          
    ##  [70] digest_0.6.37           timechange_0.3.0        R6_2.6.1               
    ##  [73] textshaping_1.0.1       colorspace_2.1-1        lpSolve_5.6.23         
    ##  [76] dichromat_2.0-0.1       biomaRt_2.58.0          mMCPcounter_1.1.0      
    ##  [79] RSQLite_2.4.3           diptest_0.77-2          tidyr_1.3.1            
    ##  [82] generics_0.1.4          data.table_1.17.8       class_7.3-23           
    ##  [85] prettyunits_1.2.0       httr_1.4.7              S4Arrays_1.2.0         
    ##  [88] pkgconfig_2.0.3         gtable_0.3.6            blob_1.2.4             
    ##  [91] XVector_0.42.0          OmnipathR_3.10.1        htmltools_0.5.8.1      
    ##  [94] carData_3.0-5           scales_1.4.0            png_0.1-8              
    ##  [97] corrplot_0.95           knitr_1.50              tzdb_0.5.0             
    ## [100] rjson_0.2.23            checkmate_2.3.3         nlme_3.1-168           
    ## [103] curl_7.0.0              proxy_0.4-27            cachem_1.1.0           
    ## [106] GlobalOptions_0.1.2     stringr_1.5.2           testit_0.13            
    ## [109] immunedeconv_2.1.0      rootSolve_1.8.2.4       KernSmooth_2.23-26     
    ## [112] parallel_4.3.3          AnnotationDbi_1.64.1    pillar_1.11.0          
    ## [115] grid_4.3.3              logger_0.4.0            vctrs_0.6.5            
    ## [118] ggpubr_0.6.1            car_3.1-3               dbplyr_2.5.1           
    ## [121] xtable_1.8-4            evaluate_1.0.5          readr_2.1.5            
    ## [124] magick_2.8.6            mvtnorm_1.3-3           cli_3.6.5              
    ## [127] locfit_1.5-9.12         compiler_4.3.3          rlang_1.1.6            
    ## [130] crayon_1.5.3            ggsignif_0.6.4          labeling_0.4.3         
    ## [133] classInt_0.4-11         mclust_6.1.1            stringi_1.8.7          
    ## [136] psych_2.5.6             BiocParallel_1.36.0     Biostrings_2.70.1      
    ## [139] limSolve_2.0.1          Matrix_1.6-5            hms_1.1.3              
    ## [142] bit64_4.6.0-1           KEGGREST_1.42.0         statmod_1.5.0          
    ## [145] igraph_2.1.4            broom_1.0.9             memoise_2.0.1          
    ## [148] bit_4.6.0               Giotto_3.3.2            readxl_1.4.5
