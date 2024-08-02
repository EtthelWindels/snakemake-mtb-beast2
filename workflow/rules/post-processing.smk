# Plot posterior parameter distributions

# BDMM results

def _get_all_bdmm_sampl025_logs(wildcards):
    files = expand(
        "results/analyses/bdmm_sampl025/{analysis}/{analysis}.log",
        analysis = [a for a in config["run"]["bdmm_sampl025"]])
    return files

rule plot_bdmm_sampl025:
    input:
        input_logs = _get_all_bdmm_sampl025_logs
    output:
        png = "results/figures/bdmm_sampl025.png",
        svg = "results/figures/bdmm_sampl025.svg"
    script:
        "../scripts/plot_bdmm_sampl025_main.R"

def _get_all_bdmm_sampl05_logs(wildcards):
    files = expand(
        "results/analyses/bdmm_sampl05/{analysis}/{analysis}.log",
        analysis = [a for a in config["run"]["bdmm_sampl05"]])
    return files

rule plot_bdmm_sampl05:
    input:
        input_logs = _get_all_bdmm_sampl05_logs
    output:
        png = "results/figures/bdmm_sampl05.png",
        svg = "results/figures/bdmm_sampl05.svg"
    script:
        "../scripts/plot_bdmm_sampl05_suppl.R"

def _get_all_bdmm_sampl1_logs(wildcards):
    files = expand(
        "results/analyses/bdmm_sampl1/{analysis}/{analysis}.log",
        analysis = [a for a in config["run"]["bdmm_sampl1"]])
    return files

rule plot_bdmm_sampl1:
    input:
        input_logs = _get_all_bdmm_sampl1_logs
    output:
        png = "results/figures/bdmm_sampl1.png",
        svg = "results/figures/bdmm_sampl1.svg"
    script:
        "../scripts/plot_bdmm_sampl1_suppl.R"

# BDMM BDSKY results

def _get_all_bdmm_bdsky_logs(wildcards):
    files = expand(
        "results/analyses/bdmm_bdsky_fixedRe/{analysis}/{analysis}.log",
        analysis = [a for a in config["run"]["bdmm_bdsky_fixedRe"]])
    return files

rule plot_bdmm_bdsky_fixedRe:
    input:
        input_logs = _get_all_bdmm_bdsky_logs
    output:
        png = "results/figures/bdmm_bdsky_fixedRe.png",
        svg = "results/figures/bdmm_bdsky_fixedRe.svg"
    script:
        "../scripts/plot_bdmm_bdsky_suppl.R" 

# BDMM inverse migration rate prior results

def _get_all_bdmm_inverse_logs(wildcards):
    files = expand(
        "results/analyses/bdmm_inverse/{analysis}/{analysis}.log",
        analysis = [a for a in config["run"]["bdmm_inverse"]])
    return files

rule plot_bdmm_inverse:
    input:
        input_logs = _get_all_bdmm_inverse_logs
    output:
        png = "results/figures/bdmm_inverse.png",
        svg = "results/figures/bdmm_inverse.svg"
    script:
        "../scripts/plot_bdmm_inv-migr-prior_suppl.R"

# BDMM clock results

def _get_all_bdmm_clock_logs(wildcards):
    files = expand(
        "results/analyses/bdmm_clockgrid/{analysis}/{analysis}.log",
        analysis = [a for a in config["run"]["bdmm_clockgrid"]])
    return files

rule plot_bdmm_clockgrid:
    input:
        input_logs = _get_all_bdmm_clock_logs
    output:
        png = "results/figures/bdmm_clockgrid.png",
        svg = "results/figures/bdmm_clockgrid.svg"
    script:
        "../scripts/plot_bdmm_clockgrid_main.R"

# BD vs. BDMM comparisons         

def _get_all_bd_sampl025_logs(wildcards):
    files = expand(
        "results/analyses/bd_sampl025/{analysis}/{analysis}.log",
        analysis = [a for a in config["run"]["bd_sampl025"]])
    return files

rule plot_bd_bdmm_comparison:
    input:
        input_logs_bd = _get_all_bd_sampl025_logs,
        input_logs_bdmm = _get_all_bdmm_sampl025_logs
    output:
        Re_png = "results/figures/bd-bdmm_comparison_Re.png",
        Re_svg = "results/figures/bd-bdmm_comparison_Re.svg",
        time_until_transm_png = "results/figures/bd-bdmm_comparison_time_until_transm.png",
        time_until_transm_svg = "results/figures/bd-bdmm_comparison_time_until_transm.svg",
        infperiod_tot_png = "results/figures/bd-bdmm_comparison_infperiod_tot.png",
        infperiod_tot_svg = "results/figures/bd-bdmm_comparison_infperiod_tot.svg",
        transrate_infperiod_png = "results/figures/bd-bdmm_comparison_transrate_infperiod.png",
        transrate_infperiod_svg = "results/figures/bd-bdmm_comparison_transrate_infperiod.svg"
    script:
        "../scripts/plot_bd-bdmm_comparison_suppl.R"