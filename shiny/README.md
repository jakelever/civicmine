## Instructions for running Shiny app

This directory contains the code for the Shiny viewer for CIViCmine that is currently deployed at [http://bionlp.bcgsc.ca/civicmine/](http://bionlp.bcgsc.ca/civicmine/).

To get it running, follow the steps below:

1. Download this code either by cloning the entire repo, or downloading the [master ZIP file](https://github.com/jakelever/civicmine/archive/refs/heads/master.zip).
2. Open the civicmine code in RStudio. 
   - Navigate to this directory (shiny/) if needed using the file viewer in the bottom-right of RStudio
   - Set the current working directory to the Shiny folder (under More -> Set As Working Directory in the file viewer)
3. Install prerequisite R packages
   - Run `install.packages(c("shiny", "data.table", "DT", "plotly", "plyr", "dplyr", "reshape2", "RColorBrewer"))` in the R terminal
4. Run the **updateCIViC.R** script which should fetch the latest release of CIViC, adapt it and do some checks
   - Open the file and then press `Ctrl + Alt + R` to run all the R commands in the file
5. Copy in the **civicmine_collated.tsv.gz** and **civicmine_sentences.tsv.gz** files from the latest release of CIViCmine available from [Zenodo](https://doi.org/10.5281/zenodo.1472826).
   - They should be in the same directory as the app.R file
6. Unzip the civicmine files (**civicmine_collated.tsv.gz** and **civicmine_sentences.tsv.gz**). 
   - The app loads the uncompressed files. You can unzip them with a tool like gunzip in the terminal or 7Zip on Windows
   - Terminal command: `gunzip civicmine_collated.tsv.gz civicmine_sentences.tsv.gz`
   - This should create files: **civicmine_collated.tsv** and **civicmine_sentences.tsv**
7. Open the **app.R** file and select the `Run App` option at the top right of the code window.
