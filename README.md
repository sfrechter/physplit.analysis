# physplit.analysis
Packaged and versioned Analysis functions for Shahar's cells

# Install

```r
# install devtools if required
if (!require("devtools")) install.packages("devtools")

# then install direct depenendencies
devtools::install_github("jefferis/gphys")
devtools::install_github("jefferislab/physplitdata")
devtools::install_github("sfrechter/physplit.analysis")
```

The only complication is that you need to have a GitHub PAT (Personal Access Token) 
setup to install the `physplitdata` repository since it is not yet shared publicly. 
This is descrbied at https://github.com/jefferislab/physplitdata.
You can also install from a local checkout if you happen to have one – something 
like this.

```r
devtools::install_git("~/dev/R/physplitdata")
```
