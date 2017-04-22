library(TwoSampleMR)
library(dplyr)

toggle_dev("test")

ao <- available_outcomes()


# Select most reliable studies for diseases and risk factors
# Exclude lipids
traits <- filter(ao, 
		category %in% c("Disease", "Risk factor"),
		(grepl("Mixed", population) | grepl("European", population)),
		mr == 1,
		sample_size > 1000,
		nsnp > 95000,
		access %in% c("immunobase_users", "public")
	) %>%
	arrange(desc(year), desc(nsnp)) %>%
	distinct(trait, .keep_all=TRUE) %>%
	filter(! id %in% c(299:302))




# get instruments

l <- list()
m <- list()
for(i in 1:nrow(traits))
{
    message(i, " of ", nrow(traits))
    l[[i]] <- try(extract_instruments(traits$id[i], p1=5e-6))
    m[[i]] <- try(extract_outcome_data(l[[i]]$SNP, traits$id))
}

save(l, m, file="extract_everything.rdata")

