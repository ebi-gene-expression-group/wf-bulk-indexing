import os
import glob

# atom: set grammar=python:

TYPES = ['annotations', 'array_designs', 'go_ens', 'interpro', 'reactome']
bioentities_directories_to_stage = set()

def get_version(source):
    if source == 'ensembl':
        return config['ens_eg_version']
    else:
        return config['wbsp_version']

def get_bioentities_directories_to_stage():
    """
    List all the directories from the main atlas bioentities that need to be staged
    to be able to run this in a per organism level.

    The web application code running these processes descides on the species to
    run based on the files it find in the BIOENTITIES path given.
    """
    global bioentities_directories_to_stage
    global TYPES
    species = config['species']
    if bioentities_directories_to_stage:
        return bioentities_directories_to_stage
    dirs=set()
    for source in ["ensembl", "wbsp"]:
        for type in TYPES:
            dir = f"{config['bioentities_source']}/archive/{type}_{get_version(source)}"
            if os.path.isdir(dir):
                dirs.add(dir)
                continue
            dir = f"{config['bioentities_source']}/archive/{type}_{source}_{get_version(source)}"
            if os.path.isdir(dir):
                dirs.add(dir)
                continue
            dir = f"{config['bioentities_source']}/archive/{type}"
            if os.path.isdir(dir):
                dirs.add(dir)

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
    if 'go_ens' in dir:
        return f"{prefix}/go"
    for type in TYPES:
        if type in dir:
            return f"{prefix}/{type}"

def get_all_staging_files():
    species = config['species']
    results = []
    source_dirs = get_bioentities_directories_to_stage()
    for sdir in source_dirs:
        dest = get_destination_dir(sdir)
        files = [os.path.basename(f) for f in glob(f"{sdir}/*{species}*")]
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
            command=f"""
                    dest={get_destination_dir(dir)}
                    mkdir -p \$dest
                    rsync -a --include '*{params.species}*' {dir}/ \$dest
                    """
            shell(command)
            print(f"{dir} staged")


rule run_bioentities_JSONL_creation:
    container: "docker://quay.io/ebigxa/atlas-index-base:1.0"
    log: "create_bioentities_jsonl.log"
    input:
        staged_files=rules.stage_files_for_species.output.staged_files
    params:
        bioentities="./",
        output_dir=config['output_dir'],
        atlas_env_file=config['atlas_env_file'],
        experiment_files="./experiment_files",
        atlas_exps=config['atlas_exps'],
        exp_design_path=config['atlas_exp_design']
    output:
        jsonl=get_jsonl_path()
    shell:
        """
        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export output_dir={params.output_dir}
        export EXPERIMENT_FILES={params.experiment_files}
        export BIOENTITIES_JSONL_PATH={params.output_dir}

        {workflow.basedir}/index-bioentities/bin/create_bioentities_jsonl.sh
        """

rule delete_species_bioentities_index:
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.0"
    params:
        atlas_env_file=config['atlas_env_file'],
        species=config['species']
    input:
        jsonl=rules.run_bioentities_JSONL_creation.output.jsonl
    output:
        deleted=touch(f"{config['species']}.index.deleted")
    shell:
        """
        source {params.atlas_env_file}
        export SPECIES={params.species}

        {workflow.basedir}/index-bioentities/bin/delete_bioentities_species.sh
        """

rule load_species_into_bioentities_index:
    container:
        "docker://quay.io/ebigxa/atlas-index-base:1.0"
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
        source {params.atlas_env_file}

        export BIOENTITIES={params.bioentities}
        export EXPERIMENT_FILES={params.experiment_files}
        export BIOENTITIES_JSONL_PATH={params.output_dir}

        {workflow.basedir}/index-bioentities/bin/index_organism_annotations.sh
        """
