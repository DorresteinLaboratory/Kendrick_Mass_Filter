# Kendrick Mass Filter
Notebooks which accompany material submitted for publication [link placed here]. *link currently non-functional until publication

**Authors:** Ricardo R. da Silva (ridasilva@ucsd.edu), Madeleine Ernst (mernst@ucsd.edu), Alan K. Jarmusch (ajarmusch@ucsd.edu) <br>
**Verion:** 1.0 (Prior to Submission of Manuscript) <br>
**Date of Last Revision:** 06-21-2018 <br>

# Description
Kendrick mass (https://en.wikipedia.org/wiki/Kendrick_mass) is the mass-to-charge (*m/z*) of each ion rescaled to the Kendrick mass scale which differs from IUPUAC mass. The scaling is performed using the nomimal mass of the unit repeat divided by the monoisotopic mass of the unit repeat. The defect between the Kendrick scaled m/z and the integer Kendrick mass value, i.e. Kendrick mass defect (KMD), is similar between homologous compounds. The Kendrick Mass Filter (KMF) is used to perform selection and removal of data centered around a user-defined KMD, with addtional restrictions on the elution time and presence of possibly multiple homologous compounds eluting together. This notebook calculates and visualizes data output obtained through the Kendrick Mass Filter for a selected dataset using user-defined parameters.

# Inputs
**Feature table:** (.csv) file with MS features in columns and samples in rows. Feature IDs are provided in the column names in the following format: "mz;RT". The first column must contain sample names.

# Output 
**Summary tables:** (.csv) files containing Kendrick mass filtered output data.
**Plots (.pdf):** Kendrick mass plot, MS1 features plot and spectra before and after applying KMF.

# Dependencies
R version 3.4.2 (2017-09-28) and packages: tidyr_0.8.0, dplyr_0.7.4, gridExtra_2.3, gtable_0.2.0, Rgraphviz_2.22.0, graph_1.56.0, BiocGenerics_0.24.0, Hmisc_4.0-3, ggplot2_2.2.1, Formula_1.2-2, survival_2.41-3, lattice_0.20-35
