# GSOC-BayesBD

To install the BayesBD package you will need:
  1.  R, version 3.4.1 was used to build BayesBD.
  2.  Rtools.
  3.  The devtools package.
  
Steps:
  1.  Open Rgui.
  2.  Install the "devtools" package if you have not already using the dropdown "Packages -> Install package(s)...".  Load the devtools      package using the command library(devtools).  
  3.  Submit the command: install_github("nasyring/GSOC-BayesBD", subdir = "BayesBD").
  
Alternative Installation:
  Thie version 1.2 is built for Windows on CRAN.  Install using the menu in Rgui.
  
Learning about the package:
  1.  Load the package with the command libary(BayesBD).  Bring up documentation using the command ?BayesBD, ?fitBinImage, ?fitContImage.
  2.  Run the examples therein.
  3.  Run BayesBDshiny() to bring up an interactive, html-based version of the package.
