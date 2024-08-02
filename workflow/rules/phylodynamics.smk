# Phylodynamics     

rule beast:
    input:
        alignment = lambda wildcards: "".join(["data/alignments","/",_get_analysis_param(wildcards, "dataset"),"/",_get_analysis_param(wildcards, "dataset"),"_beast.fasta"]),
        xml = lambda wildcards: _get_analysis_param(wildcards, "xml"),
    output:
        trace = "results/analyses/{method}/{analysis}/chains/{analysis}.{chain}.log",
        trees = "results/analyses/{method}/{analysis}/chains/{analysis}.{chain}.trees",
    params:
        beast_command = lambda wildcards: _get_analysis_param(wildcards, "command"), 
        action = lambda wildcards: _get_analysis_param(wildcards, "action"),
        user_check = lambda wildcards: _get_analysis_param(wildcards, "user_check"),
        xml_params = lambda wildcards: str(_get_analysis_param(wildcards, "xml_params")).replace(":", "=").replace(
            "{", "\"").replace("}", "\"").replace(" ", "").replace("'", ""),
        folder_name = "results/analyses/{method}/{analysis}/chains",
        file_name = "{analysis}.{chain}",
        state_file = "{analysis}.{chain}.state",
    log:
        "logs/{analysis}_{chain}.txt"
    benchmark:
        "benchmarks/{analysis}_{chain}.benchmark.txt"
    threads:
        lambda wildcards: _get_analysis_param(wildcards, "threads"),
    resources:
        runtime = lambda wildcards: _get_analysis_param(wildcards, "time"),
        mem_per_cpu = lambda wildcards: _get_analysis_param(wildcards, "mem_per_cpu")
    shell:
        """
        mkdir -p {params.folder_name}/running
        
        # for resuming, make sure the log and tree files are in the running folder

        {params.beast_command} \
            -D aligned_fasta={input.alignment} \
            -D {params.xml_params} \
            -D file_name={params.folder_name}/running/{params.file_name} \
            -seed {wildcards.chain} \
            -statefile "{params.folder_name}/{params.state_file}" \
            -{params.action} {input.xml} 2>&1 | tee -a {log}
        
        # log and tree files need to be written to a different path than snakemake output so snakemake does not delete them if job fails; move files once job is finished:
        [ -f {params.folder_name}/running/{params.file_name}.log ] && mv {params.folder_name}/running/{params.file_name}.log {output.trace}
        [ -f {params.folder_name}/running/{params.file_name}.trees ] && mv {params.folder_name}/running/{params.file_name}.trees {output.trees} 
        """

def _get_chains(wildcards):
    files = expand(
        "results/analyses/{{method}}/{{analysis}}/chains/{{analysis}}.{chain}.{{output}}",
        chain = range(1, _get_analysis_param(wildcards, "chains") + 1))
    return files

rule combine_chains:
    message: 
        """
        Combine chain files: {input.chain_files} with LogCombiner.
        """
    input:
        chain_files = _get_chains  
    output:
        combined_chain = "results/analyses/{method}/{analysis}/{analysis}.{output}",
    params:
        burnin =  lambda wildcards: _get_analysis_param(wildcards, "burnin"),
        resample =  lambda wildcards: _get_analysis_param(wildcards, "resample"),
        input_command = lambda wildcards, input: " -log ".join(input) 
    resources:
        mem_per_cpu = 4096
    shell:
        """
        logcombiner -log {params.input_command} -o {output.combined_chain} -b {params.burnin}  -resample {params.resample} 2>&1 | tee -a {log}
        """