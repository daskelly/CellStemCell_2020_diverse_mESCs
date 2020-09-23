## Obtaining targets of hotspots

In the manuscript we identified several eQTL and caQTL hotspots.
We can obtain the targets of each hotspot by loading some
information from Supplementary Table 2.

Table S2 consists of the "Genome Coordinates of Hotspot QTL" (worksheet 1)
and a "List of Mediation Results for QTL Peaks That May Be Regulated by Each 
Hotspot" (worksheet 2). For the mediation analysis we considered each gene
located within or < 5Mb from the hotspot boundaries to be a candidate
mediator and used each in a mediation analysis with each target.


```r
library(tidyverse)
library(readxl)
```

Switch to a directory where you've downloaded Table S2.
Let's get the targets of the Chr 10 eQTL hotspot.

```r
hotspot <- "Chr 10, 61 Mb"
dat <- read_excel("Table_S2.xlsx", sheet=2, 
    col_types=c('text', 'text', 'text', 'text', 'numeric', 'numeric', 'text', 'text', 'numeric'))
targets <- filter(dat, hotspot_id == hotspot) %>% select(target_id, target_symbol) %>% distinct()
```

Note that, depending on the hotspot, 
the exact number of hotspot targets might differ slightly
from numbers discussed in some sections of the manuscript. Initial
reports of hotspots were based on mapping conducted independently
in the gene expression or chromatin accessibility datasets,
as appropriate. 
In contrast, Table S2 presents results from mediation analyses,
which were performed using the subset of samples present in both the 
gene expression and chromatin accessibility datasets.

