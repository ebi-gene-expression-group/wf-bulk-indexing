import os
import glob

# atom: set grammar=python:

TYPES = ['annotations', 'array_designs', 'go', 'interpro', 'reactome']
bioentities_directories_to_stage = set()

print(f"ENS version: {config['ens_version']}")
print(f"ENS GN version: {config['eg_version']}")

def get_version(source):
    if source == 'ensembl' or source == "ens":
        return f"{config['ens_version']}_{config['eg_version']}"
    else:
        return config['wbsp_version']

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
    prefix=f"{config['bioentities_source']}/archive"
    for type in TYPES:
        if type == "array_designs":
            dir = f"{prefix}/{type}_{get_version('ens')}_{config['wbsp_version']}"
            if os.path.isdir(dir):
                print(f"{dir} exists")
                dirs.add(dir)
                continue
        for source in ["ensembl", "wbps", "ens"]:
            versionSep = "_"
            if source == "ens":
                versionSep = ""
            dir = f"{prefix}/{type}_{get_version(source)}"
            if os.path.isdir(dir):
                print(f"{dir} exists")
                dirs.add(dir)
                continue
            dir = f"{prefix}/{type}_{source}{versionSep}{get_version(source)}"
            if os.path.isdir(dir):
                print(f"{dir} exists")
                dirs.add(dir)
                continue
            dir = f"{prefix}/{type}"
            if os.path.isdir(dir):
                print(f"{dir} exists")
                dirs.add(dir)
                continue
            dir = f"{prefix}/{source}_{get_version(source)}"
            if os.path.isdir(dir):
                print(f"{dir} exists")
                dirs.add(dir)
                continue

    bioentities_directories_to_stage = dirs
    return bioentities_directories_to_stage

def get_destination_dir(dir):
    """
    Provides the directories where the files will be staged to
    for the bioentities properties, based on the source directory.

    Assumes that the worklfow runs in the temporary bioentities
    created for the species.
    """
    prefix='bioentity_properties'
    for type in TYPES:
        if type in dir:
            return f"{prefix}/{type}"
    for source in ["ensembl", "wbps"]:
        if source in dir:
            return f"{prefix}/{source}"

def get_all_staging_files():
    species = config['species']
    results = []
    source_dirs = get_bioentities_directories_to_stage()
    for sdir in source_dirs:
        dest = get_destination_dir(sdir)
        print(f"Looking at {sdir} with destination {dest}")
        files = [os.path.basename(f) for f in glob.glob(f"{sdir}/{species}*")]
        if dest.endswith("go"):
            files = [os.path.basename(f) for f in glob.glob(f"{sdir}/*")]
        results.extend([f"{dest}/{f}" for f in files])

    print(results)
    return results

def get_jsonl_path():
    return f"{config['output_dir']}/{config['species']}.jsonl"



rule stage_files_for_species:
    log: "staging.log"
    input:
        directories=get_bioentities_directories_to_stage()
    output:
        staged_files=get_all_staging_files()
    params:
        species=config['species']
    run:
        for dir in input.directories:
            dest = get_destination_dir(dir)
            if dest.endswith("go"):
                call = f"rsync -a {dir}/* {dest}"
            else:
                # some directories which are not "go" will not have anything for our species
                if not glob.glob(f"{dir}/{params.species}*"):
                    print(f"Skipping {dir} for {params.species}")
                    continue
                call = f"rsync -a {dir}/{params.species}* {dest}"

            print(f"Calling {call}")
            command = f"""
                      mkdir -p {dest}
                      {call}
                      """
            shell(command)
            print(f"{dir} staged")


rule run_bioentities_JSONL_creation:
    container: "docker://quay.io/ebigxa/atlas-index-base:1.0"
    log: "create_bioentities_jsonl.log"
    input:
        #staged="staged.done"
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
        jsonl=get_jsonl_path()
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

        mkdir -p {params.experiment_files}
        rm -f {params.experiment_files}/magetab
        rm -f {params.experiment_files}/expdesign
        ln -s {params.atlas_exps} {params.experiment_files}/magetab
        ln -s {params.exp_design_path} {params.experiment_files}/expdesign
        ln -s {params.web_app_context}/species-properties.json {params.experiment_files}/species-properties.json
        ln -s {params.web_app_context}/release-metadata.json {params.experiment_files}/release-metadata.json

        if [ -f /bin/micromamba ]; then
            eval "$(/bin/micromamba shell hook -s bash)"
            micromamba activate "$ENV_NAME"
        fi

        {workflow.basedir}/index-bioentities/bin/create_bioentities_jsonl.sh
        """

rule delete_species_bioentities_index:
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.0"
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

        if [ -f /bin/micromamba ]; then
            eval "$(/bin/micromamba shell hook -s bash)"
            micromamba activate "$ENV_NAME"
        fi

        {workflow.basedir}/index-bioentities/bin/delete_bioentities_species.sh
        """

rule load_species_into_bioentities_index:
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.0"
    log: "load_species_into_bioentities_index.log"
    params:
        bioentities="./",
        output_dir=config['output_dir'],
        atlas_env_file=config['atlas_env_file'],
        experiment_files="./experiment_files",
        atlas_exps=config['atlas_exps'],
        exp_design_path=config['atlas_exp_design']
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

        {workflow.basedir}/index-bioentities/bin/index_organism_annotations.sh
        """
