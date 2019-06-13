#!/bin/sh

# Loop over years.
   YEARS=("1992:2000" "2001:2009" "2010:2018" "1992:2018")
# Loop over regression indices.
   REG=("none" "NAO" "PDO")

# Do loop.
   for Y_IDX in {0..3}; do
      for R_IDX in {0..2}; do
         sed -e 's/<anal_years>/'${YEARS[${Y_IDX}]}'/'  \
             -e 's/<flag_reg>/'${REG[${R_IDX}]}'/'      \
         CCI_trend_alldata_mpi_template.R > CCI_trend_mpi_${Y_IDX}_${R_IDX}.R

# Run R script.
         mpirun -np 4 Rscript CCI_trend_mpi_${Y_IDX}_${R_IDX}.R &> mpi_out_${Y_IDX}_${R_IDX}
         mv CCI_trend_mpi_${Y_IDX}_${R_IDX}.R run_files
      done
   done

