
NBCONVERT = jupyter-nbconvert --ExecutePreprocessor.timeout=600 --to notebook --execute

.SECONDARY:

.PHONY: stats all

# Macro/function for getting data path for a given experiment's sample dataset
sampling_data_path = data/$(1)/$(1)_sampling_results.csv
scheme_base = $(subst MARC_,,$(1))
get_schemes = $(call scheme_base,$(1))_3 \
			  $(call scheme_base,$(1))_4 \
			  ARG MBN

STATS = MARC_main_stats_processed.p MARC_gCCN_stats_processed.p

all: oneone oneone_aerosol cdf stats_tables sobol_lhs

stats: stats_tables

stats_tables: $(STATS)
	$(NBCONVERT) stats_tables.ipynb
MARC_main_stats_processed.p: data/MARC_main/MARC_main_sampling_results.csv
	python compute_stats.py $< $(call get_schemes,MARC_main)
MARC_gCCN_stats_processed.p: data/MARC_gCCN/MARC_gCCN_sampling_results.csv
	python compute_stats.py $< $(call get_schemes,MARC_gCCN)

# Set of oneone sampling plots
oneone:
	$(NBCONVERT) oneone_plots.ipynb

# Set of oneone sampling plots, for MARC aerosols
oneone_aerosol:
	$(NBCONVERT) oneone_MARC_aerosol_plots.ipynb

# CDF comparison plots
cdf:
	$(NBCONVERT) cdf_plots.ipynb

# Sobol indices comparison for LHS sampling
sobol_lhs:
	$(NBCONVERT) sobol_decomp.ipynb

clean_nbs:
	rm *.nbconvert.ipynb
clean:
	rm oneone cdf stats_tables $(STATS) *.nbconvert.ipynb
