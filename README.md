# GNPS
Camil Gosmanov

github.com/camilgosmanov

Description:
Find mirror matches from mass spectrometry data uploaded to GNPS using R

1. Download cytoscape files from the GNPS job status page ("Advanced Views - External Visualization" -> "Direct Cytoscape Preview/Download")
2. Open in Cytoscape and export the Edge table and Node table to downloads as .csv
3. Make the features table as .csv with M/Z_RT as column 1
4. Rename the Edge table, Node table, and features table files to "Edge.csv", "Nodes.csv", and "Features.csv", or just manually change the path in the program
5. Set the url to the GNPS job status page
6. Set the working directory and output file names, then run the program (GNPS Matches + Mirror Plots.R)
7. When the browser closes, the program is finished and the .pptx and .csv file will be output to the working directory

Notes:
- This program uses firefox to navigate through the GNPS library
- If you are using Windows, you will need to install Rtools to make the presentation
- If you get an httr error, you may need to download the latest version of either Java or Java JDK
- If you get a library error, you may need to download the latest version of either R or R Studio 
