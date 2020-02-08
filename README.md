# CAAs_1
Chromosome arm aneuploidies (CAAs) in cancer

This README file pertains to:

Chromosome arm aneuploidies shape tumour evolution and drug response.
Shukla A, Nguyen THM, Moka SB, Ellis JJ, Grady JP, Oey H, Cristino AS, Khanna KK, Kroese DP, Krause L, Dray E, Fink JL, 
Duijf PHG.
Nature Communications 11: 449 (2020).
https://www.nature.com/articles/s41467-020-14286-0

This pan-cancer project involved analyses of chromosome arm aneuploidies (CAAs) in cancer development. The study of CAAs 
revealed novel aspects of tumour evolution, metastasis, patient survival outcome, and chemotherapeutic drug response for a 
broad range of cancers.

This repository contains data files and Python and R scripts used in the study.

SYSTEM REQUIREMENTS
- Mac OS (OS X El Capitan or later) or Microsoft Windows (Windows 10 or later)
- R (version 3.6 or later)
- Python (version 3 or later)

INSTALLATION GUIDE
- R software: Refer to: 
"https://cran.r-project.org/bin/windows/base/" (Mac OS X), 
"https://cran.r-project.org/bin/windows/base/" (Windows).
Install time: 1-2 minutes.
- R packages: Refer to .R R scripts provided in the "scripts" folder.
- Python 3 software: Refer to: 
"https://www.python.org/downloads/macos" (Mac OS X), 
"https://www.python.org/downloads/windows/" (Windows).
Install time: 1-2 minutes.
- Python modules and packages: Refer to .py Python scripts provided in the "scripts" folder.

DEMO
- Install Python/R software, modules and packages, as indicated above and in provided scripts.
- Import modules/install packages/load libraries as indicated in scripts.
- Demo to call CAAs from source data:
  * Use input file and Python script in the "CAAs_1/demo" folder. Segmented CN data are in the subfolder "CAAs_1/demo/cnv_data". (Unzip file "MSK-IMPACT_CNA_Segments.seg.zip" first.)
  * Type at the shell prompt: "$ python3 MSKCM-IMPACT_CAL-SNCAs_cnv_analysis.py".
  * The python script uses the .seg file for processing. Note that the input file for the segmented CN data inside the "cnv_data" folder needs to have the extension ".seg" and that the segmented CN data are tab-delimited and in UTF-8 format.
- Run other scripts (in "CAAs_1/scripts" folder) on data (in "CAAs_1/data files" folder or the "Source Data" file).
- Expected output: See Supplemental Tables accompanying the paper.

INSTRUCTIONS FOR USE
- Install software as described under INSTALLATION GUIDE above.
- Run software as per example described under DEMO above.
