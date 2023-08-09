if [ -f "../output/demo_results_overview.pdf" ]; then
	rm ../output/demo_results_overview.pdf
fi

# Rscript slim_sims_pipeline.R bottle 1000 1000 50 250 10000 10900 10950 11000
# Rscript slim_sims_pipeline.R small 1000 250 250 250 10000 10300 10600 11000
# Rscript slim_sims_pipeline.R decline 1000 1000 250 50 10000 10900 10950 11000
Rscript slim_sims_pipeline.R large 1000 2000 2000 2000 10000 10300 10600 11000

Rscript demo_results_quick_check.R
