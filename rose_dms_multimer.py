import os

configfile: "config.yaml"
WTP     = config["wt_pdb"]
WTPREF  = os.path.splitext(os.path.basename(WTP))[0]
OUTDIR  = config["outdir"]
XML     = config["xml"]
ROSETTA = config["rosetta_bin"]
LD_PATH = config["ld_path"]

POSITIONS = config.get("positions", [68, 79, 81, 194, 254, 256, 470, 471, 516, 517, 641])
AAS = ["ALA","ARG","ASN","ASP","CYS","GLU","GLN","GLY","HIS","ILE",
       "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL"]
WT_AAS = {int(k): v for k, v in config["wt_aas"].items()}


rule all:
    input:
        os.path.join(OUTDIR, config["ddg_csv"]),
        os.path.join(OUTDIR, config["rmsd_csv"])


rule rosetta_mutscan:
    input:
        pdb = WTP,
        xml = XML
    output:
        pdb = os.path.join(OUTDIR, f"{WTPREF}_{{pos}}_{{aa}}_0001.pdb"),
        sc  = os.path.join(OUTDIR, f"{WTPREF}_{{pos}}_{{aa}}_score.sc")
    params:
        prefix = os.path.join(OUTDIR, f"{WTPREF}_{{pos}}_{{aa}}_"),
        ld     = LD_PATH
    shell:
        """
        export LD_LIBRARY_PATH={params.ld}:$LD_LIBRARY_PATH
        {ROSETTA} \
            -s {input.pdb} \
            -parser:protocol {input.xml} \
            -parser:script_vars focused_res={wildcards.pos} target_aa={wildcards.aa} \
            -out:prefix {params.prefix} \
            -nstruct 1
        """


rule calc_ddg:
    input:
        expand(os.path.join(OUTDIR, f"{WTPREF}_{{pos}}_{{aa}}_score.sc"),
               pos=POSITIONS, aa=AAS)
    output:
        os.path.join(OUTDIR, config["ddg_csv"])
    params:
        score_dir = OUTDIR
    shell:
        "python3 calc_ddg.py {params.score_dir} {output}"


rule calc_rmsd:
    input:
        pdbs = expand(os.path.join(OUTDIR, f"{WTPREF}_{{pos}}_{{aa}}_0001.pdb"),
                      pos=POSITIONS, aa=AAS)
    output:
        os.path.join(OUTDIR, config["rmsd_csv"])
    params:
        score_dir = OUTDIR
    shell:
        "python3 calc_rmsd.py {params.score_dir} {output}"
