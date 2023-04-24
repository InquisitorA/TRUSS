# Compute_Truss
Compute truss of a graph by distributing it among different processors and running them in parallel

OpenMPI facilitates communication between processors. Massive graphs are taken as input and read in chunks separately by different processes, first support is calculated for each edge and the graph is filtered and truss is computed. Additionaly OpenMP is used to speed up computations within a single processor by running several threads in parallel.
