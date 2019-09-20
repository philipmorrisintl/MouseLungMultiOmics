# Multi-omics systems toxicology study of mouse lung tissue assessing the biological effects of aerosols from two heat-not-burn tobacco products and cigarette smoke

This repository contains the R analysis code and R data objects for
the analysis of lung multi-omics data reported in Titz et
al. (submitted).

## Content

* **SCRIPTS/P15038_APOE_P2_MultiOmicsManuscript.Rmd** : Rmd file with
  the analysis code
* **DATA/** : Folder with data files for each omics modality
* **DATA/EXTERNAL/** : Folder with external data files supporting the
  analysis (see below how to obtain external files)
* **INFO/** : Folder with additional annotation files

## Installation

### R

We have used R version 3.5.1, a more recent version of R should work
but hasn't been tested.

### Install CRAN packages:

```{r}
 req_packages <- c("knitr",
                   "gridExtra",
                   "RColorBrewer",
                   "ggplot2",
                   "egg",
                   "reshape2",
                   "xlsx",
                   "readxl",
                   "openxlsx",
                   "tools",
                   "plotrix",
                   "gdata",
                   "plyr",
                   "dplyr",
                   "stringi",
                   "ggbeeswarm",
                   "visNetwork",
                   "devtools")
 req_packages <- req_packages[!req_packages %in% rownames(installed.packages())]
 if (length(req_packages) > 0) {
    install.packages(req_packages)
 }
```

### Install Bioconductor packages:

```{r}
 if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
 BiocManager::install()
    
 req_packages <- c("limma",
                   "MOFA",
                   "mixOmics")
 
 req_packages <- req_packages[!req_packages %in% rownames(installed.packages())]
 if (length(req_packages) > 0) {
    BiocManager::install(req_packages)
 }
```

Consult vignette of MOFA package for further details on the
installation, including the setup of the python environment.
 
### Install github packages:

```{r}
 if (!"PCSF" %in% rownames(installed.packages())) {
    BiocManager::install("topGO")
    devtools::install_github("IOR-Bioinformatics/PCSF", repos=BiocManager::repositories(),
                             dependencies=TRUE, type="source", force=TRUE)
 }
    
 if (!"NPA" %in% rownames(installed.packages())) {
    devtools::install_github("philipmorrisintl/NPAModels",
                             dependencies=TRUE, 
                             type="source")
    devtools::install_github("philipmorrisintl/NPA",
                             dependencies=TRUE, 
                             type="source")
 }
 
 #available from Bioconductor (see above)
 #devtools::install_github("bioFAM/MOFA", build_opts = c("--no-resave-data"))
 
```

### Python

Any Python3 version should be ok, we used 3.6.4.
Please note Python < 2.7 is not supported.

* Create a Python virtualenv

You can create it anywhere you have access.
```bash
$ python3 -m venv .
```

Activate it and install the necessary `mofapy`.
```bash
$ source bin/activate

$ pip install -U pip setuptools
(...)
Successfully installed pip-19.1.1 setuptools-41.0.1

$ pip install mofapy
(...)
Successfully installed argparse-1.4.0 h5py-2.9.0 joblib-0.13.2 mofapy-1.2 numpy-1.16.4 pandas-0.24.2 python-dateutil-2.8.0 pytz-2019.1 scikit-learn-0.21.2 scipy-1.3.0 six-1.12.0 sklearn-0.0
```
Please note the version numbers are just an example.

At any time you can deactivate the venv with the following command.
```bash
$ deactivate
```
### Download this code/data package from github

Create a new folder for the project which could be the same as the
Python's virtualenv, but this is not required, download and unzip
the repository (example, if done from R environment):

```{r}
#set destination folder
project_folder = "path/to/project/folder"

#create folder
dir.create(project_folder)
setwd(project_folder)

#download from github
download.file(url = "https://github.com/philipmorrisintl/MouseLungMultiOmics/archive/1.0.1/1.0.1.zip", 
              destfile = "Lung_MultiOmics.zip")
              
#unzip
unzip(zipfile = "Lung_MultiOmics.zip")

#list content
list.files()

```

### Obtain gene-set collections & network files
* mSigDB
    + Check license terms before download
    + Download gene-set collection files (gmt format) from
      http://software.broadinstitute.org/gsea/downloads.jsp
    + c2.all.v6.2.symbols.gmt and h.v6.2.symbols.gmt are required
    + Save both files in DATA/EXTERNAL folder of project

* StringDB
    + Check license terms before download
    + Download functional interaction network files from
      https://version-10-5.string-db.org/cgi/download.pl?species_text=Mus+musculus
    + 10090.protein.aliases.v10.5.txt and
      10090.protein.links.detailed.v10.5.txt are required (version
      10.5, unzip required)
    + Save both txt files in DATA/EXTERNAL folder
    
* miRTarBase
    + Check license terms before download
    + Download
      http://mirtarbase.mbc.nctu.edu.tw/cache/download/7.0/mmu_MTI.xls
    + We downloaded the version of 10 July 2018 (version 7.0)
    + Save in DATA/EXTERNAL folder with name
      mirTarBase_Mm_10July2018.xls (or adjust file name in script)
      
* Reactome
    + Check license terms before download
    + Download und unzip https://reactome.org/download/current/ReactomePathways.gmt.zip
    + We downloaded the version of 23 Apr 2018
    + Save gmt file as 'ReactomePathways_23Apr2018.gmt'
      in DATA/EXTERNAL folder (or adjust file name in script)
    
* KEGG (optional)
    + Special license required see: https://www.kegg.jp/kegg/download/
    + Download und unzip kegg/genes/organisms/mmu/mmu_link.tar.gz
    + Download und unzip kegg/genes/organisms/mmu/T01002.kff.gz
    + Download und unzip kegg/ligand/reaction.tar.gz
    + Save reaction_ko.list, reaction_mapformula.lst, and T01002.kff
      in DATA/EXTERNAL folder
      
## Run analysis
* Activate the previously created Python virtual environment
* Start R and change directory to the main project folder
  (the folder with this README file).
* Adjust run options in script as necessary (force_rerun_mofa,
  force_rerun_pcsf, force_recreate_network). Note that
  a KEGG license is required to obtain the corresponding files,
  creating the integrated network, and running the PCSF analysis.
* Execute the analysis and render the PDF file with the following command:

```{r}
    rmarkdown::render("SCRIPTS/P15038_APOE_P2_MultiOmicsManuscript.Rmd", 
                      output_dir = "REPORT",
                      intermediates_dir = "REPORT",
                      clean = FALSE)

```

The generated files, figures, and the rendered PDF report are located
in the **REPORT/** folder.

## License

We are releasing the analysis code (SCRIPTS/ folder) under the following license:

> Copyright 2019 Philip Morris Products SA
> 
> This program is free software; you can redistribute it and/or
> modify it under the terms of the GNU General Public License
> as published by the Free Software Foundation; either version 2
> of the License, or (at your option) any later version.
> 
> This program is distributed in the hope that it will be useful,
> but WITHOUT ANY WARRANTY; without even the implied warranty of
> MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
> GNU General Public License for more details.
> 
> You should have received a copy of the GNU General Public License
> along with this program; if not, write to the Free Software
> Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
  
  
Also, see Notices.txt for the licenses of the used packages/libraries.
  
  
For the shared data (DATA/ & INFO/ folders) :

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

## References
Titz et al. Multi-omics systems toxicology study of mouse lung tissue assessing the biological effects of 
aerosols from two heat-not-burn tobacco products and cigarette smoke. *submitted*

  
## Related Repositories
- NPA (Network Perturbation Amplitude): <https://github.com/philipmorrisintl/NPA>
- NPAModels (R package and data for NPA models): <https://github.com/philipmorrisintl/NPAModels>

## Contact
Bjoern Titz (<DL.RSupport@pmi.com>)


