#!/bin/sh

# Season.
   MONTHS="JFMOND"
# Loop over years.
   YEARS=("1992-2000" "2001-2009" "2010-2018" "1992-2018")
# Loop over regression indices.
   REG=("none" "NAO" "PDO")

# Do loop.
   for Y_IDX in {0..2}; do
      for R_IDX in {0..2}; do
         sed -e 's/<anal_years>/'${YEARS[${Y_IDX}]}'/'  \
             -e 's/<flag_reg>/'${REG[${R_IDX}]}'/'      \
             -e 's/<lab_months>/'${MONTHS}'/'      \
         ggplot_CCI_trend_NEA_template.R > ggplot_CCI_NEA_${Y_IDX}_${R_IDX}.R

# Run R script.
         Rscript ggplot_CCI_NEA_${Y_IDX}_${R_IDX}.R &> ggplot_out_${Y_IDX}_${R_IDX}
         mv ggplot_CCI_NEA_${Y_IDX}_${R_IDX}.R run_files
      done
   done

