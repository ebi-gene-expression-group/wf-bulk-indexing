#!/usr/bin/env bash

# This script takes the marker_gene files (normally available on
# staging) for a gxa experiment and does the following steps:
# - Transforms those files into a long (melted) table of experiment id,
#   cell/run id, and expression. This takes care of avoiding large chunks of
#   data being kept in memory for long, at the expense of writing to disk (too)
#   often.
# - Deletes data from gxa_marker_gene table on postgres if exists 
# - Loads data into the gxa_marker_gene table.
# 
set -e

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}" )" &> /dev/null && pwd )

postgres_scripts_dir="${SCRIPT_DIR}/../postgres_routines"

checkDatabaseConnection() {
  psql $1 -c '\q' > /dev/null 2>&1 || { echo "PostgreSQL is not ready"; }
}

dbConnection=${dbConnection:-$1}
EXP_ID=${EXP_ID:-$2}
metrices="tpms fpkms" # add proteomics

# Check that necessary environment variables are defined.
[ -z ${dbConnection+x} ] && echo "Env var dbConnection for the database connection needs to be defined. This includes the database name." && exit 1
[ -z ${EXP_ID+x} ] && echo "Env var EXP_ID for the id/accession of the experiment needs to be defined." && exit 1

# Check that files are in place.
for metric in $metrices;
do
    marker_genes_path=$ATLAS_PROD/analysis/baseline/*/experiments/${EXP_ID}/${EXP_ID}-${metric}-markers.tsv
    
    ls $marker_genes_path > /dev/null 2>&1 || {
        echo "No matching file found for $marker_genes_path"
        exit 1
    }

done

# Check that database connection is valid
checkDatabaseConnection $dbConnection
    
# Deletes existing marker genes from gxa_marker_gene table 
sed "s/<EXP-ACCESSION>/$EXP_ID/" $postgres_scripts_dir/01-delete_existing_marker_gene.sql.template | \
psql -v ON_ERROR_STOP=1 "$dbConnection"

# Load gene marker table
    
for metric in $metrices;
do
    marker_genes_path=$ATLAS_PROD/analysis/baseline/*/experiments/${EXP_ID}/${EXP_ID}-${metric}-markers.tsv

    echo $marker_genes_path

    marker_genes_path=$(ls $marker_genes_path)

    echo $marker_genes_path

    # removes header
    tail -n +2 ${marker_genes_path} > ${marker_genes_path}.no_header.tsv
    
    sed "s|<PATH-TO-DATA>|${marker_genes_path}.no_header.tsv|" $postgres_scripts_dir/02-load_gene_marker_table.sql.template | \
    psql -v ON_ERROR_STOP=1 $dbConnection

    rm -rf ${marker_genes_path}.no_header.tsv
done
