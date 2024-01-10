import os
import glob
import re
from textwrap import dedent
from pathlib import Path

# atom: set grammar=python:

# Note about Solr and db authentication:
# The scripts used in this workflow expect SOLR_USER and SOLR_PASS environment variables
# set to authenticate. This is currently expected to be set in the ATLAS_ENV_FILE sourced 
# at every task. Since the same var names are used in the env file and on the scripts,
# there is no need for them to be made explicit in each rule. The same happens with db auth.

TYPES = ['annotations', 'array_designs', 'go', 'interpro', 'reactome', 'mirbase']
bioentities_directories_to_stage = set()
staging_files = set()

micromamba_env = dedent("""
                 if [ -f /bin/micromamba ]; then
                     eval "$(/bin/micromamba shell hook -s bash)"
                     micromamba activate "$ENV_NAME"
                 fi
                 """)

get_missing_accessions = dedent("""
                         # needs $prefix (a directory) defined
                         export failed_accessions_output=$prefix"/failed_accessions.txt"
                         done_accessions=$prefix"/done_accessions.txt"
                         if [ -f $done_accessions ]; then
                             #find lines only in input_accessions (that haven't been done)
                             comm -23 <( sort -u $input_accessions ) <( sort -u $done_accessions ) > $prefix"/missing_accessions.txt"
                             input_accessions=$prefix"/missing_accessions.txt"
                         fi
                         """)


def get_bioentities_directories_to_stage():
    """
    List all the directories from the main atlas bioentities that need to be staged
    to be able to run this in a per organism level.

    The web application code running these processes descides on the species to
    run based on the files it find in the BIOENTITIES path given.

    This is the structure of directories that we aim to match
    annotations_ensembl_104_51*  annotations_wbps_15      go_ens104_51*
    array_designs_104_51_15*     ensembl_104_51*          reactome_ens104_51*  wbps_15*
    """
    global bioentities_directories_to_stage
    global TYPES
    species = config['species']
    print(f"Lenght of set {len(bioentities_directories_to_stage)}")
    if bioentities_directories_to_stage:
        return bioentities_directories_to_stage
    dirs=set()
    prefix=f"{config['bioentities_source']}"
    for type in TYPES:
        dir = f"{prefix}/{type}"
        if os.path.isdir(dir):
            print(f"{dir} exists")
            dirs.add(dir)

    bioentities_directories_to_stage = dirs
    return bioentities_directories_to_stage


def get_destination_dir(dir):
    prefix='bioentity_properties'
    return f"{prefix}/{os.path.basename(dir)}"

def get_all_staging_files():
    global staging_files
    if staging_files:
        return staging_files
    species = config['species']
    results = []
    source_dirs = get_bioentities_directories_to_stage()
    for sdir in source_dirs:
        dest = get_destination_dir(sdir)
        print(f"Looking at {sdir} with destination {dest}")
        if dest.endswith("go") or dest.endswith('interpro'):
            files = [os.path.basename(f) for f in glob.glob(f"{sdir}/*.tsv")]
        else:
            files = [os.path.basename(f) for f in glob.glob(f"{sdir}/{species}.*.tsv")]
            # We have cases where the files are buried one directory below :-(
            files.extend([os.path.basename(f) for f in glob.glob(f"{sdir}/*/{species}.*.tsv")])
        results.extend([f"{dest}/{f}" for f in files])

    staging_files.update(results)
    return staging_files

wbps_annotations_re = re.compile(r"annotations.*wbps")
ens_annotations_re = re.compile(r"annotations.*ens")
array_designs_re = re.compile(r"array_designs.*\.(A-[A-Z]{4}-[0-9]+)\.tsv")

def get_jsonl_label(input):
    global wbps_annotations_re
    global ens_annotations_re
    global array_designs_re

    if "reactome" in input:
        return "reactome"
    if "mirbase" in input:
        return "mature_mirna"
    if wbps_annotations_re.search(input):
        return "wbpsgene"
    if ens_annotations_re.search(input):
        return "ensgene"
    m_array = array_designs_re.search(input)
    if m_array:
        return m_array.group(1)


def get_jsonl_paths():
    inputs_for_jsonl = get_all_staging_files()
    jsonls = set()
    for input in inputs_for_jsonl:
        json_label = get_jsonl_label(input)
        if json_label:
            jsonls.add(f"{config['output_dir']}/{config['species']}.{json_label}.jsonl")
    print(f"Number of JSONLs expected: {len(jsonls)}")
    return jsonls

def species_for_db(species):
    """
    Atlas database and web app use
    species as "Homo sapiens" instead of "homo_sapiens".
    """
    return species.replace("_", " ").capitalize()

def get_mem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    """
    mem_avail = [ 2, 2, 4, 8, 16, 64, 128, 256 ]
    return mem_avail[attempt-1] * 1000

# get_mem_mb_coexp
def get_himem_mb(wildcards, attempt):
    """
    To adjust resources in the rules
    attemps = reiterations + 1
    Max number attemps = 8
    Beyond 300 GBs it needs a different queue on our setup, but snakemake
    currently only supports using a single queue.
    """
    mem_avail = [ 8, 8, 16, 32, 64, 128, 256, 290 ]
    return mem_avail[attempt-1] * 1000

def get_coexp_mem_mb(wildcards, attempt):
    mem_avail = [8, 32, 64, 128, 160, 192, 208, 240]
    return mem_avail[attempt-1] * 1000

def get_sync_cpus(wildcards):
    return 4

def get_exp_design_sync_destination():
    # The destination for a sync as expected by rsync: user@host:<path>
    if 'exp_update_sync_dest' in config:
        return config['exp_update_sync_dest']


def aggregate_accessions_update_experiment(wildcards):
    checkpoint_output = checkpoints.divide_accessions_into_chunks.get(**wildcards).output[0]
    return expand("update_experiment_designs/{chunk}/exp_designs_updates.txt",
        chunk=glob_wildcards("accessions_{chunk}").chunk)

def aggregate_accessions_sync_experiment_designs(wildcards):
    checkpoint_output = checkpoints.divide_accessions_into_chunks.get(**wildcards).output[0]
    return expand("update_experiment_designs/{chunk}/sync_experiment_designs.done",
        chunk=glob_wildcards("accessions_{chunk}").chunk)

def aggregate_baseline_accessions_update_coexpression(wildcards):
    checkpoint_output = checkpoints.divide_accessions_into_chunks.get(**wildcards).output[0]
    return expand("update_coexpressions/{chunk}/update_coexpressions.txt",
        chunk=glob_wildcards("baseline_accessions_{chunk}").chunk)

def aggregate_accessions_load_bulk_analytics_index(wildcards):
    checkpoint_output = checkpoints.divide_accessions_into_chunks.get(**wildcards).output[0]
    return expand("load_bulk_analytics_index/{chunk}/analytics_index_loaded.txt",
        chunk=glob_wildcards("accessions_{chunk}").chunk)

localrules: aggregate_update_experiment, aggregate_update_coexpression, aggregate_load_bulk_analytics_index


rule get_accessions_for_species:
    """
    This step brings all accessions for a desired species, it leaves all the accessions (including baseline) 
    in `species_accessions.txt` and all the baseline accessions in `species_baseline_accessions.txt`. This is
    useful for the update setup (after doing an E! Update).

    When you want to use this same pipeline for specific accessions loading you will need to:
    
    1.- Obtain all the accessions that need reloading
    2.- Obtain for each accession its species (and check that those species are in bioentities)
    3.- Group the accessions per species (each species will be a separate run).
    4.- For each species, create the files `species_accessions.txt` and `species_baseline_accessions.txt` 
        only with the desired accesisons.

    Step 4 should ensure that this specific rule doesn't run, and that only the desired accessions are loaded.
    """
    log: "get_accessions_for_species.log"
    input:
        exclusion_file=config['exclude_accessions']
    params:
        species=species_for_db(config['species']),
        atlas_env_file=config['atlas_env_file']
    output:
        accessions=temp("species_accessions.txt"),
        baseline_accessions=temp("species_baseline_accessions.txt")
    conda: "envs/postgres.yaml"
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {params.atlas_env_file}

        psql -c "COPY (SELECT accession FROM experiment WHERE species LIKE '{params.species}%' ORDER BY load_date) TO STDOUT WITH NULL AS ''" \
             -v ON_ERROR_STOP=1 $dbConnection > {output.accessions}
        psql -c "COPY (SELECT accession FROM experiment WHERE species LIKE '{params.species}%' AND type LIKE '%BASELINE%' ORDER BY load_date) TO STDOUT WITH NULL AS ''" \
             -v ON_ERROR_STOP=1 $dbConnection > {output.baseline_accessions}

        echo {input.exclusion_file}
        # Loop through each word in exclude.txt
        while IFS= read -r accession && [ -n "$accession" ]; do
            # Remove the word from {output.accessions} using grep
            grep -v "$accession" {output.accessions} > temp && mv temp {output.accessions}
            grep -v "$accession" {output.baseline_accessions} > temp && mv temp {output.baseline_accessions}
        done < {input.exclusion_file}
        """

checkpoint divide_accessions_into_chunks:
    log: "divide_accessions_into_chunks.log"
    params:
        lines_per_split=50,
        lines_per_split_baseline=5
    input:
        accessions=rules.get_accessions_for_species.output.accessions,
        baseline_accessions=rules.get_accessions_for_species.output.baseline_accessions
    output:
        done=touch("divide_accessions_into_chunks.done"),
    shell:
        """
        # This will generate accessions_01, accessions_02, etc and the same for baseline_accessions_
        split -l {params.lines_per_split} -d {input.accessions} accessions_
        split -l {params.lines_per_split_baseline} -d {input.baseline_accessions} baseline_accessions_
        """

rule stage_files_for_species:
    """
    This rule brings files that are needed for the bioentities indexing of the species. These files
    will be usually at $ATLAS_PROD/bioentities (they are the result of the E! Update process) into local
    directories of this pipeline execution.
    
    This needs to happen during a general (en-masse) update (after E! Update), but not 
    during the loading of specific experiments, as no new species are dealt with at that point.
    """
    log: "staging.log"
    input:
        directories=get_bioentities_directories_to_stage()
    output:
        staged_files=get_all_staging_files()
    params:
        species=config['species']
    run:
        rsync_options = '-a --delete'
        for dir in input.directories:
            dest = get_destination_dir(dir)
            # even if there is no content, the web app context will expect the dir.
            Path(dest).mkdir(parents=True, exist_ok=True)
            if dest.endswith("go") or dest.endswith("interpro"):
                call = f"rsync {rsync_options} --include=*.tsv  --exclude=* {dir}/* {dest}"
            elif not glob.glob(f"{dir}/{params.species}.*.tsv") and not glob.glob(f"{dir}/*/{params.species}.*.tsv"):
                print(f"Skipping {dir} for {params.species}")
                continue
            elif dest.endswith('annotations') or dest.endswith('array_designs'):
                # some directories which are not "go" will not have anything for our species
                call = f"rsync {rsync_options} {dir}/**/{params.species}.*.tsv {dest}"
            elif dest.endswith('reactome') or dest.endswith('mirbase'):
                call = f"rsync {rsync_options} {dir}/{params.species}.*.tsv {dest}"


            print(f"Calling {call}")
            command = f"""
                      exec &> "{log}"
                      mkdir -p {dest}
                      {call}
                      """
            shell(command)
            print(f"{dir} staged")


rule prepare_directories_and_links:
    """
    This rule does some minor reformatting of directories after the files have been copied over
    in the stage_files_for_species rule. This again is needed everytime that we do an E! Update.
    """
    log: "prepare_directories_and_links.log"
    input:
        staged_files=rules.stage_files_for_species.output.staged_files
    params:
        bioentities="./",
        output_dir=config['output_dir'],
        atlas_env_file=config['atlas_env_file'],
        experiment_files="./experiment_files",
        atlas_exps=config['atlas_exps'],
        web_app_context=config['web_app_context'],
        exp_design_path=config['atlas_exp_design']
    output:
        dirs_prepared=touch("dirs_prepared")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {params.atlas_env_file}

        mkdir -p {params.experiment_files}
        mkdir -p {params.output_dir}
        rm -f {params.experiment_files}/magetab
        rm -f {params.experiment_files}/expdesign
        ln -sf {params.atlas_exps} {params.experiment_files}/magetab
        ln -sf {params.exp_design_path} {params.experiment_files}/expdesign
        ln -sf {params.web_app_context}/species-properties.json {params.experiment_files}/species-properties.json
        ln -sf {params.web_app_context}/release-metadata.json {params.experiment_files}/release-metadata.json
        touch {output.dirs_prepared}
        """

rule update_experiment_designs:
    """
    This step updates the experiment design files that the web application tomcat process expects in
    some local location of the machine running tomcat. These files are based on the Atlas metadata
    SDRF files. On the tomcat server these will usually end up in /srv/gxa/data/gxa/expdesign

    Happens per accession.

    Later on, rule sync_experiment_designs does an rsync to the specified tomcat machine (based on config).
    """
    container: "docker://quay.io/ebigxa/atlas-index-base:1.5"
    log: "update_experiment_designs/{chunk}/update_experiment_designs.log"
    resources:
        mem_mb=get_mem_mb
    params:
        bioentities="./",
        output_dir="update_experiment_designs/{chunk}",
        atlas_env_file=config['atlas_env_file'],
        experiment_files="./experiment_files",
        atlas_exps=config['atlas_exps'],
        web_app_context=config['web_app_context'],
        exp_design_path=config['atlas_exp_design']
    input:
        accessions="accessions_{chunk}",
        dirs_prepared=rules.prepare_directories_and_links.output.dirs_prepared
    output:
        done=touch("update_experiment_designs/{chunk}/exp_designs_updates.txt")
    shell:
        """
        prefix={params.output_dir}
        mkdir -p $prefix
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export output_dir={params.output_dir}
        export EXPERIMENT_FILES={params.experiment_files}
        export server_port=8081 #fake
        export PS1="\\u@\\h:\\w\\$ "


        input_accessions={input.accessions}

        {micromamba_env}

        {get_missing_accessions}

        # only run if we have accessions to run (>0)
        if [ $(wc -l $input_accessions | awk '{{ print $1 }}') -gt 0 ]; then
            export ACCESSIONS=$(cat $input_accessions | tr '\\n' ',' | sed 's/,$//' )
            set +e
            {workflow.basedir}/index-gxa/bin/update_experiment_designs_cli.sh
            status=$?
            touch $failed_accessions_output # in case there was no failure.
            # check if the process failed but the failed file is empty
            # if this is the case, the process failed before starting and we
            # shouldn't update the done accessions
            if ( [[ -s $failed_accessions_output ]] && ((status != 0)) ) || (( status == 0 )); then
                # accessions done now, input accessions that are not in failed -- append to all done accessions
                comm -23 <( sort -u $input_accessions ) <( sort -u $failed_accessions_output ) >> $done_accessions
            fi
            set -e
            exit $status
        else
            echo "No more accessions to be done on this retry."
        fi
        """


rule sync_experiment_designs:
    """
    Usually experiment designs are stored in the local disk mounted by the web application server. Consider that when setting
    config['exp_update_sync_dest']. This step doesn't get executed if this config is not set.
    """
    log: "update_experiment_designs/{chunk}/sync_experiment_designs.log"
    threads: get_sync_cpus
    input:
        source_accs="accessions_{chunk}",
        update_exp_designs_done="update_experiment_designs/{chunk}/exp_designs_updates.txt"
    params:
        destination=get_exp_design_sync_destination(),
        exp_design_prefix="experiment_files/expdesign"
    output:
        done=touch("update_experiment_designs/{chunk}/sync_experiment_designs.done")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set

        cat {input.source_accs} | parallel -j {threads} --joblog {log} rsync -qlrtvz --no-owner {params.exp_design_prefix}/ExpDesign-{{}}.tsv {params.destination}/
        """


rule update_coexpressions:
    """
    Coexpressions are currently stored in the database. This process uses the web app cli to update those. No element from this process is stored in the file system.

    Happens per accession.
    """
    container: "docker://quay.io/ebigxa/atlas-index-base:1.5"
    log: "update_coexpressions/{chunk}/update_coexpressions.log"
    resources:
        mem_mb=get_coexp_mem_mb
    params:
        bioentities="./",
        output_dir="update_coexpressions/{chunk}",
        atlas_env_file=config['atlas_env_file'],
        experiment_files="experiment_files",
        atlas_exps=config['atlas_exps'],
        web_app_context=config['web_app_context'],
        exp_design_path=config['atlas_exp_design']
    input:
        baseline_accessions="baseline_accessions_{chunk}",
        dirs_prepared=rules.prepare_directories_and_links.output.dirs_prepared
    output:
        done=touch("update_coexpressions/{chunk}/update_coexpressions.txt")
    shell:
        """
        prefix={params.output_dir}
        mkdir -p $prefix
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export output_dir={params.output_dir}
        export EXPERIMENT_FILES={params.experiment_files}
        export server_port=8081 #fake

        input_accessions={input.baseline_accessions}

        {micromamba_env}

        {get_missing_accessions}

         # only run if we have accessions to run (>0)
        if [ $(wc -l $input_accessions | awk '{{ print $1 }}') -gt 0 ]; then
            export ACCESSIONS=$(cat $input_accessions | tr '\\n' ',' | sed 's/,$//' )
            set +e
            {workflow.basedir}/index-gxa/bin/update_coexpressions_cli.sh
            status=$?
            touch $failed_accessions_output # in case there was no failure.
            # check if the process failed but the failed file is empty
            # if this is the case, the process failed before starting and we
            # shouldn't update the done accessions
            if ( [[ -s $failed_accessions_output ]] && ((status != 0)) ) || (( status == 0 )); then
                # accessions done now, input accessions that are not in failed -- append to all done accessions
                comm -23 <( sort -u $input_accessions ) <( sort -u $failed_accessions_output ) >> $done_accessions
            fi
            set -e
            exit $status
        else
            echo "No more accessions to be done on this retry."
        fi
        """


rule aggregate_update_experiment:
    input: aggregate_accessions_update_experiment
    output: "exp_designs_updates.done"
    shell:
        """
        touch {output}
        """

rule aggregate_sync_experiment_designs:
    input: aggregate_accessions_sync_experiment_designs
    output: "sync_exp_designs.done"
    shell:
        """
        touch {output}
        """

rule aggregate_update_coexpression:
    input: aggregate_baseline_accessions_update_coexpression
    output: "coexpression_updates.done"
    shell:
        """
        touch {output}
        """

rule run_bioentities_JSONL_creation:
    """
    This rule takes the output of the E! Update, which was copied over from $ATLAS_PROD/bioentities in
    the staging and prepare directories rules and generates JSONL files with the content that needs to be filled
    into the bioentities index.

    This should only happen when you are running a whole species update, but shouldn't be done
    when we are loading new experiments or reloading existing ones. The logic to avoid this run
    for specific accession loading is not currently implemented.

    TODO: avoid this to run in the per experiment loading scenario.
    """
    container: "docker://quay.io/ebigxa/atlas-index-base:1.5"
    log: "create_bioentities_jsonl.log"
    input:
        staged_files=rules.stage_files_for_species.output.staged_files,
        dirs_prepared=rules.prepare_directories_and_links.output.dirs_prepared
    params:
        bioentities="./",
        output_dir=config['output_dir'],
        atlas_env_file=config['atlas_env_file'],
        experiment_files="./experiment_files",
        atlas_exps=config['atlas_exps'],
        web_app_context=config['web_app_context'],
        exp_design_path=config['atlas_exp_design']
    resources:
        mem_mb=16000
    output:
        jsonl=get_jsonl_paths()
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export output_dir={params.output_dir}
        export EXPERIMENT_FILES={params.experiment_files}
        export BIOENTITIES_JSONL_PATH={params.output_dir}
        export server_port=8081 #fake

        {micromamba_env}

        {workflow.basedir}/index-bioentities/bin/create_bioentities_jsonl.sh
        """

rule delete_species_bioentities_index:
    """
    This step deletes the current contents of bioentities index for the specified species. We want this
    for E! Update downstream loading, but not for specific accessions loading. The logic to avoid that for the
    specific accessions run is not currently implemented.

    This won't run if the JSONLs in run_bioentities_JSONL creation rule are not created.
    """
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.5"
    log: "delete_species_bioentities_index.log"
    params:
        atlas_env_file=config['atlas_env_file'],
        species=config['species']
    input:
        jsonl=rules.run_bioentities_JSONL_creation.output.jsonl
    output:
        deleted=touch(f"{config['species']}.index.deleted")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {params.atlas_env_file}
        export SPECIES={params.species}

        {micromamba_env}

        {workflow.basedir}/index-bioentities/bin/delete_bioentities_species.sh
        """

rule load_species_into_bioentities_index:
    """
    This step loads the current contents of JSONL files into bioentities index for the specified species. We want this
    for E! Update downstream loading, but not for specific accessions loading. The logic to avoid that for the
    specific accessions run is not currently implemented.

    This won't run if the JSONL and previous deletion don't happen.
    """
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.5"
    log: "load_species_into_bioentities_index.log"
    params:
        bioentities="./",
        output_dir=config['output_dir'],
        atlas_env_file=config['atlas_env_file'],
        experiment_files="./experiment_files",
        atlas_exps=config['atlas_exps'],
        exp_design_path=config['atlas_exp_design'],
        species=config['species']
    input:
        jsonl=rules.run_bioentities_JSONL_creation.output.jsonl,
        deleted_confirmation=rules.delete_species_bioentities_index.output.deleted
    output:
        loaded=touch(f"{config['species']}.index.loaded")
    shell:
        """
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"
        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export EXPERIMENT_FILES={params.experiment_files}
        export BIOENTITIES_JSONL_PATH={params.output_dir}
        export SPECIES={params.species}
        export server_port=8081 #fake

        {micromamba_env}

        # for the json_loader script which is called in turn by index_organism_annotations.
        export PATH={workflow.basedir}/index-bioentities/bin:$PATH

        {workflow.basedir}/index-bioentities/bin/index_organism_annotations.sh
        """

######################### analytics #########################

rule analytics_bioentities_mapping:
    """
    This rule creates a binary file with mappings from specific accessions' genes to bioentities
    elements (as available in $ATLAS_PROD/bioentities). As such, this step doesn't require the JSONL
    or the data loaded into bioentities solr index, but will use the files that where staged in the 
    rules stage_files_for_species and prepare_directories_and_links.

    This mapping file is later used to fill in the analytics index, which is accession specific.

    This needs to happen both in general updating loading (post E! Update) and for accession specific loading.
    """
    log: "analytics_bioentities_mapping/{chunk}/analytics_mapping.log"
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.5"
    input:
        # This could optionally be either that file or a file given with specific accessions to redo.
        # or maybe the accessions broken in chunks.
        accessions="accessions_{chunk}",
        exp_design="update_experiment_designs/{chunk}/exp_designs_updates.txt",
        index_loaded=rules.load_species_into_bioentities_index.output.loaded
    params:
        bioentities="./",
        output_dir="analytics_bioentities_mapping/{chunk}",
        atlas_env_file=config['atlas_env_file'],
        experiment_files="experiment_files",
        species=config['species']
    output:
        mappings_file="analytics_bioentities_mapping/{chunk}/"+f"{config['species']}.map.bin"
    resources:
        mem_mb=get_himem_mb
    shell:
        """
        prefix={params.output_dir}
        mkdir -p $prefix
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export EXPERIMENT_FILES={params.experiment_files}
        export output_dir={params.output_dir}
        export SPECIES={params.species}
        export server_port=8081 #fake

        # needed to trigger an error code exit for mappings
        export failed_accessions_output=$prefix"/failed_accessions.txt"

        input_accessions={input.accessions}

        {micromamba_env}

        # we don't break this into done and failed as we need a single mapping file per chunk
        # Note that part of the data from the mapping comes from the Solr server, so as an input the fact that bioentites has been loaded is needed.

        export ACCESSIONS=$(cat $input_accessions | tr '\\n' ',' | sed 's/,$//' )
        {workflow.basedir}/index-bioentities/bin/create_bioentities_property_map.sh
        """

rule create_analytics_jsonl_files:
    """
    This rule create the JSONL files, relying on the binary mapping files, for each accession to be loaded
    later into the bulk-analytics index in Solr.

    This needs to happen both in general updating loading (post E! Update) and for accession specific loading.
    """
    log: "analytics_jsonl_files/{chunk}/analytics_jsonl_files.log"
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.5"
    input:
        # This could optionally be either that file or a file given with specific accessions to redo.
        # or maybe the accessions broken in chunks.
        accessions="accessions_{chunk}",
        mappings_file="analytics_bioentities_mapping/{chunk}/"+f"{config['species']}.map.bin"
    resources:
        mem_mb=get_himem_mb
    params:
        bioentities="./",
        output_dir="analytics_jsonl_files/{chunk}",
        mappings_directory="analytics_bioentities_mapping/{chunk}/",
        atlas_env_file=config['atlas_env_file'],
        experiment_files="experiment_files",
        species=config['species']
    output:
        created=touch("analytics_jsonl_files/{chunk}/analytics_jsonl_files.txt")
    shell:
        """
        prefix={params.output_dir}
        mkdir -p $prefix
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export EXPERIMENT_FILES={params.experiment_files}
        export output_dir={params.output_dir}
        export SPECIES={params.species}
        export server_port=8081 #fake
        export BIN_MAP={params.mappings_directory}

        input_accessions={input.accessions}

        {micromamba_env}

        {get_missing_accessions}

        mkdir -p {params.output_dir}

        # only run if we have accessions to run (>0)
        if [ $(wc -l $input_accessions | awk '{{ print $1 }}') -gt 0 ]; then
            export ACCESSIONS=$(cat $input_accessions | tr '\\n' ',' | sed 's/,$//' )
            set +e
            {workflow.basedir}/index-gxa/bin/generate_analytics_JSONL_files.sh
            status=$?
            touch $failed_accessions_output # in case there was no failure.
            # check if the process failed but the failed file is empty
            # if this is the case, the process failed before starting and we
            # shouldn't update the done accessions
            if ( [[ -s $failed_accessions_output ]] && ((status != 0)) ) || (( status == 0 )); then
                # accessions done now, input accessions that are not in failed -- append to all done accessions
                comm -23 <( sort -u $input_accessions ) <( sort -u $failed_accessions_output ) >> $done_accessions
            fi
            set -e
            exit $status
        else
            echo "No more accessions to be done on this retry."
        fi
        """

rule load_bulk_analytics_index:
    """
    This rule loads all the analytics JSONL, created in the rule create_analytics_jsonl_files,
    into the bulk-analytics-index.

    This needs to happen both in general updating loading (post E! Update) and for accession specific loading.
    """
    log: "load_bulk_analytics_index/{chunk}/load_bulk_analytics_index.log"
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.5"
    input:
        jsonl_created=rules.create_analytics_jsonl_files.output.created,
        accessions="accessions_{chunk}"
    params:
        bioentities="./",
        output_dir="load_bulk_analytics_index/{chunk}",
        analytics_jsonl_dir="analytics_jsonl_files/{chunk}",
        atlas_env_file=config['atlas_env_file'],
        experiment_files="experiment_files",
        species=config['species']
    output:
        loaded=touch("load_bulk_analytics_index/{chunk}/analytics_index_loaded.txt")
    shell:
        """
        prefix={params.output_dir}
        mkdir -p $prefix
        set -e # snakemake on the cluster doesn't stop on error when --keep-going is set
        exec &> "{log}"

        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export EXPERIMENT_FILES={params.experiment_files}
        export SPECIES={params.species}
        export server_port=8081 #fake

        input_accessions={input.accessions}

        {micromamba_env}

        {get_missing_accessions}

        export PATH={workflow.basedir}/index-gxa/bin:$PATH

        # only run if we have accessions to run (>0)
        if [ $(wc -l $input_accessions | awk '{{ print $1 }}') -gt 0 ]; then
            export ACCESSIONS=$(cat $input_accessions | tr '\\n' ',' | sed 's/,$//' )
            export delete_existing=true
            export analytics_jsonl_dir={params.analytics_jsonl_dir}
            # generates id clashes
            # {workflow.basedir}/index-gxa/bin/gxa-index-set-autocreate.sh
            set +e
            {workflow.basedir}/index-gxa/bin/load_analytics_files_in_Solr.sh
            status=$?
            touch $failed_accessions_output # in case there was no failure.
            # check if the process failed but the failed file is empty
            # if this is the case, the process failed before starting and we
            # shouldn't update the done accessions
            if ( [[ -s $failed_accessions_output ]] && ((status != 0)) ) || (( status == 0 )); then
                # accessions done now, input accessions that are not in failed -- append to all done accessions
                comm -23 <( sort -u $input_accessions ) <( sort -u $failed_accessions_output ) >> $done_accessions
            fi
            set -e
            exit $status
        else
            echo "No more accessions to be done on this retry."
        fi
        """

rule aggregate_load_bulk_analytics_index:
    input: aggregate_accessions_load_bulk_analytics_index
    output: "load_bulk_analytics_index.done"
    shell:
        """
        touch {output}
        """
