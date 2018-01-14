# Genome relations browser - test.

Explore the origin of admixed varieties following genome wide classification.

Visit the [app](https://ideogram.herokuapp.com/)

## Context

This application represents a first attempt at a more comprehensive and intuitive way of exploring relations among
individuals along their genomes.

The results presented follow a genome crawl that classifies varieties locally according to their proximity
to pre-defined reference populations. This study was conducted on rice populations (*Oryza sativa*), and members of
the *Aus*, *Indica* and *Japonica* groups, as well as a considerable number of *circum*Basmati and general admixed accessions.

The result of this crawl is the classification, along their genome, of each rice variety into *pure* or
*intermediate* classes. Pure classes - in yellow, blue, red and black - represent unambiguous classifications
into one of the reference populations, Aus, Japonica and Indica, or to an outlier class, respectively.
Unambiguous classificaitons were allowed for when assignment probabilities relative to two or three classes
did not differ by more than a factor of two (pair-wise comparisons).

While the global result of this genome crawl is already informative, allowing one to move past the
traditional attribution of genome proportions and to identify which regions were likely inherited from whom
(in populational terms at least), it is limited as well. Parent populations are structured, and one would
like to be able to identify which subpopulations are in fact closer to the ancestral introgressants.

By selecting the region and reference population to focus on, this application attempts to provide an
intuitive and fast way to do just that.

More information is provided on the application itself.

## App features
Range slider focuses analysis on specified region, updates ideogram to show only that region.

Cluster observations (graph): Principal component analysis of clusters extracted. Colors derive from a simple kmeans clustering at K == 9.

Likelihood density (graph): density plot of mean likelihood of each accessions for clusters of a given color.

Accessions - loadings (graph): Loadings plot of the above PCA. Reveals overall relationships between accessions given the clusters selected for this analysis. color-code: Blue= Japonica; yellow= circumAus; red= Indica; purple= Admix; green= cBasmati.

Selected accessions (table): Information on either all 948 CORE accessions included (select Full-colors in "Chose a color" dropdown), or those above a certain threshold when focusing on a particular cluster (select color "Chose a color", modulate threshold likelihood using slider).

## Use
**Online:**

Visit the [Application](https://ideogram.herokuapp.com/) hosted on Heroku servers.

