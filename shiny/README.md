## Instructions for running Shiny app

This directory contains the code for the Shiny viewer for CIViCmine that is currently deployed at [http://bionlp.bcgsc.ca/civicmine/](http://bionlp.bcgsc.ca/civicmine/).

To get it running, follow the steps below:

1. Download this code (either by cloning the entire repo, or downloading the [master ZIP file](https://github.com/jakelever/civicmine/archive/refs/heads/master.zip).
2. Open the shiny directory in RStudio
  - Set the current working directory to be the Shiny directory if needed
3. Install prerequisite R packages
  - Run `install.packages(c("shiny", "data.table", "DT", "R.utils", "plotly", "plyr", "dplyr", "reshape2", "RColorBrewer"))` in the R terminal
4. Run the **updateCIViC.R** script which should fetch the latest release of CIViC, adapt it and do some checks
5. Copy in the **civicmine_collated.tsv.gz** and **civicmine_sentences.tsv.gz** files from the latest release of CIViCmine available from [Zenodo](https://doi.org/10.5281/zenodo.1472826).
6. Open the **app.R** file and select the `Run App` option at the top right of the code window.