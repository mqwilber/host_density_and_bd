The file `candidate_failed_invasion_sites.csv' gives the unique lakes with R. muscosa/sierrae in the Sierra Nevada that potentially experienced a failed Bd invasion. The criteria for a candidate failed invasion is

(Bd was detected in the population at year t) AND 
(Bd was not detected in the population at year t + 1) AND 
(At least 1 swab was taken from the population at year t + 1)

This is a purposefully broad criteria and will also include lakes that experienced epidemic fade-out.

The columns in the data file are

`site_id': An integer. Unique lake identification number

`in_roland': A boolean value (TRUE or FALSE). If TRUE, the given lake was identified to have experienced a potential failed invasion in the full dataset provided by Roland on Nov. 29, 2020. If FALSE, no failed invasion was identified at this lake, but one was identified using data provided by Cherie (see below).   These data span years 2004-2019.  These data only included the first Bd swab replicate for a given swab ID.  There are 5 lakes that Cherie's data identified as having failed invasions that were not identified by Roland's data. These differences are because Cherie's data includes replicate Bd swabs for a given swab ID and in these 5 lakes Bd loads < 1 were identified for a few years (potentially false positives).

`in_cherie': A boolean value (TRUE or FALSE). If TRUE, the given lake was identified to have experienced a potential failed invasion in the dataset provided by Cherie. These data span years 2004-2014.  These data include multiple Bd swab replicates for a given swab ID.  Other the 5 lakes mentioned above, the lakes identified as having candidate failed invasions by Cherie's data are a subset of those in Roland's data (which is what we could expect).

`loss_year`: A character string.  Gives the yearly transition(s) where a lake went from Bd positive to Bd negative.  Multiple  1 -> 0 transitions may have occurred in a single lake. These are separated by commas e.g. (2009 -> 2010, 2011 -> 2012).

`num_swabs': A character string. The number of Bd swabs taken on frogs at a lake in the year where Bd was not detected (the second year in loss_year).  Multiple numbers are given if multiple candidate failed invasions where observed.

`area': An integer. The area of the lake

`jurisdiction': A character string. The administrative location of the lake within the Sierra Nevada.

`failed_invasion_ranking':  An integer. FOR TOM AND ROLAND TO FILL IN. Gives Tom and Roland's observation on how likely it is that the candidate lake actually experienced a failed invasion.
    1 - Unlikely that a failed invasion occurred at this site
    2 - Could have experienced a failed invasion
    3 - Likely experienced a failed invasion

