# GNPS
Description:
Find mirror matches from mass spectrometry data uploaded to GNPS using R (must have firefox installed)

1. Download cytoscape files from the GNPS job status page ("Advanced Views - External Visualization" -> "Direct Cytoscape Preview/Download")
2. Open in Cytoscape and export the Edge table and Node table to downloads as .csv
3. Make the features table as .csv with M/Z_RT and fold change columns
4. Rename the Edge table, Node table, and features table files to "Edge.csv", "Nodes.csv", and "Matches.csv", or just manually change the path in the program
5. Set the url to the GNPS job status page
6. Set the working directory and output file names, then run the program
7. When program is running in firefox, keep the window open in front to minimize loading time errors
8. When the browser closes, the program is finished and the .pptx and .csv file will be output to downloads
