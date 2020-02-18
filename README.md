==NEM Pipeline==

To generate simulated data with false positive and false negative rates as specified, run methods, and plot results.  Results will be written into the directory specified.

=== Installation===

```
devtools::install_github("bitmask/labnetmet")
devtools::install_github("bitmask/NEMpipeline")
```


Optional dependencies:
```
lem: install.packages("https://www.mimuw.edu.pl/~szczurek/lem/lem_1.0.tar.gz",repos=NULL, type="source")
bnem: devtools::install_github("MartinFXP/B-NEM")
pcnem: devtools::install_github("cbg-ethz/pcNEM")
fgnem: https://sysbio.soe.ucsc.edu/projects/fgnem/
```


```
library(NEMpipeline)
set.seed(42)
alpha <- 0.15
beta <- 0.05
output_dir <- "~/projects/NEMpipelineoutput"
run_simulated_pipeline(alpha, beta, output_dir)
```
