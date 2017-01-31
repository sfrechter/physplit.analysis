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
This is described at https://github.com/jefferislab/physplitdata.

## Problems
If you run into an error along the lines of:

```
Downloading GitHub repo jefferislab/physplitdata@master
from URL https://api.github.com/repos/jefferislab/physplitdata/zipball/master
Error in stop(github_error(request)) : Not Found (404)
```
at the last step, then you can replaced the last step with. :
```r
devtools::install_github("sfrechter/physplit.analysis", dependencies=FALSE)
```
See e.g. https://github.com/hadley/devtools/issues/1381 for the origin of this error.
## Hacking physplit.analysis
You can clone this repo and use RStudio to work with `physplitdata.Rproj` to 
work on the functions in this package.

You can also install from a local checkout if you happen to have one – something 
like this.

```r
devtools::install_git("~/dev/R/physplitdata")
```
