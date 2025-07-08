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

error_exit() {
  echo "[ERROR] $1" >&2
  exit 1
}


SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]:-$0}")" &>/dev/null && pwd)"
POSTGRES_SCRIPTS_DIR="${SCRIPT_DIR}/../postgres_routines"
METRICS=("tpms" "fpkms")  # Add more metrics as needed

EXP_ID=${EXP_ID:-$1}

# Check that necessary environment variables are defined.
[ -z ${dbConnection+x} ] && error_exit  "Env var dbConnection for the database connection needs to be defined."
[ -z ${EXP_ID+x} ] && error_exit  "Env var EXP_ID for the id/accession of the experiment needs to be defined."


# Check that database connection is valid
check_db_connection() {
  psql "$dbConnection" -c '\q' &>/dev/null || error_exit "PostgreSQL is not ready or connection failed."
  echo "Database connection established using $dbConnection."
}
    
# Deletes existing marker genes from gxa_marker_gene table 
delete_old_data() {
  local sql_file="${POSTGRES_SCRIPTS_DIR}/01-delete_existing_marker_gene.sql.template"
  [[ ! -f "$sql_file" ]] && error_exit "SQL template not found: $sql_file"
  sed "s/<EXP-ACCESSION>/$EXP_ID/" "$sql_file" | psql -v ON_ERROR_STOP=1 "$dbConnection"
}

find_marker_file() {
  local metric="$1"
  local path_pattern="${ATLAS_EXPS}/${EXP_ID}/${EXP_ID}-${metric}-markers.tsv"
  local file
  file=$(ls $path_pattern 2>/dev/null | head -n1)
  echo "$file"
}

# Load gene marker table
load_marker_data() {
  local metric="$1"
  local marker_file
  marker_file=$(find_marker_file "$metric")
  echo "Processing file: $marker_file"

  local no_header_file="${marker_file}.no_header.tsv"
  tail -n +2 "$marker_file" > "$no_header_file"

  local sql_file="${POSTGRES_SCRIPTS_DIR}/02-load_gene_marker_table.sql.template"
  [[ ! -f "$sql_file" ]] && error_exit "SQL template not found: $sql_file"
  sed "s|<PATH-TO-DATA>|$no_header_file|" "$sql_file" | psql -v ON_ERROR_STOP=1 "$dbConnection"

  rm -f "$no_header_file"
}


check_db_connection
delete_old_data

found_any=0
for metric in "${METRICS[@]}";
do
    marker_file=$(find_marker_file "$metric")
    if [[ -n "$marker_file" ]]; then
      load_marker_data "$metric"
      found_any=1
    else
        echo "No marker file found for metric: ${metric} in ${EXP_ID}"
    fi
done

if [[ $found_any -eq 0 ]]; then
    error_exit "No marker files found for metrics ${METRICS[@]} in ${EXP_ID}."
fi
