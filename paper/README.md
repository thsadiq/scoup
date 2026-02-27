# Folder Information
This folder contains files necessary to facilitate submission of an article that describes the `scoup` package to the Journal of Open Source Software.

The [`data`](data) sub-folder contains the image used in the manuscript and the necessary files and outputs that produced the image. These include the following.

- [molecules](data/molecules): a sub-folder that contains all the protein sequences generated using the `scoup` package. The files are named in terms of the values of the variance of the synonymous (`000 = 0.00`, `010=0.01`) and the non-synonymous (`001=0.001`, `005=0.005`, `009=0.009`,  `030=0.030`,  `070=0.070`,  `100=0.100`, `500=0.500`,  `900=0.900`) selection coefficients. For example, the file named `s000n500.txt` was simulated with a non-synonymous coefficient variance equal to 0.50 and a synonymous coefficient variance equal to 0.00. The evolutionary tree that was used during the simulation is also saved in this sub-folder as [`simtree.txt`](data/molecules/simtree.txt).

- [procfiles](data/procfiles): another sub-folder that contains the control `.ctl` files that were used to implement PAML to obtain the maximum likelihood estimates of the ratio of non-synonymous to synonymous substitution rates for the simulated data in the [molecules](data/molecules) sub-folder. The `.ctl` files are named in the same way as the data they were used to process.

- [s000dnds](data/s000dnds.csv): an Excel file that contains the analytical estimates of the ratio of non-synonymous to synonymous substitution rates for the simulated data where the synonymous coefficient variance was set to 0.00. The values were obtained as outputs from the `scoup` simulation.

- [s010dnds](data/s010dnds.csv): an Excel file that contains the analytical estimates of the ratio of non-synonymous to synonymous substitution rates for the simulated data where the synonymous coefficient variance was set to 0.01. The values were obtained as outputs from the `scoup` simulation.

- [s000omega](data/s000omega.csv): an Excel file that contains the maximum likelihood estimates of the ratio of non-synonymous to synonymous substitution rates for the simulated data where the synonymous coefficient variance was set to 0.00. The values were obtained as outputs from PAML.

- [s010omega](data/s010omega.csv): an Excel file that contains the maximum likelihood estimates of the ratio of non-synonymous to synonymous substitution rates for the simulated data where the synonymous coefficient variance was set to 0.01. The values were obtained as outputs from PAML.

- [simulate](data/simulate.R): contains the `R` code that was used to execute the simulation that generated the data in the [molecules](data/molecules) sub-folder and the [s000dnds](data/s000dnds.csv) \& [s010dnds](data/s010dnds.csv) files.

- [summarise](data/summarise.R): contains the `R` code that was used to process the analyses output and generate the results presented in the manuscript.

