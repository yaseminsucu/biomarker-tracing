import numpy as np
import pandas as pd
import os, logging, shutil, subprocess, re, time
import scipy.stats as stats
from typing import List


GENE_ID_SYMBOLS = "/mnt/vm-shared-storage/Gene_ID_symbols_files/gene_id_symbol_df.tsv"
GENE_ID_HGNC = "/mnt/vm-shared-storage/Gene_ID_symbols_files/gene_id_symbol_hgnc.tsv"



def load_prot_data(base_path: str, disease: str, atlas: pd.DataFrame):
    gene_names_other = pd.read_csv(GENE_ID_SYMBOLS, sep="\t").rename(columns={"gene_ids": "gene", "gene_symbols": "gene_name"})
    gene_names_hgnc = pd.read_csv(GENE_ID_HGNC, sep="\t")
    gene_names_hgnc = gene_names_hgnc.drop(columns=["Status", "Approved_name", "HGNC_ID", "NCBI_gene_ID", "UCSC_gene_ID"]).rename(columns={"Approved_symbol": "gene_name", "Ensembl_gene_ID": "gene"})
    gene_names = pd.concat([gene_names_other, gene_names_hgnc], axis=0).drop_duplicates()

    mappings = {
        "BAP18": "ENSG00000258315", "CERT": "ENSG00000113163", "GPR15L": "ENSG00000188373",
        "KIR2DL2": "ENSG00000273661", "HLA": "ENSG00000204592", "MENT": "ENSG00000143443",
        "LEG1": "ENSG00000184530", "LILRA3": "ENSG00000278046", "NTproBNP": "ENSG00000120937",
        "HLA-DRA": "ENSG00000204287", "PALM2": "ENSG00000157654", "SARG": "ENSG00000182795"
    }

    prot_df = pd.read_csv(f"{os.path.join(base_path, disease)}.csv")
    prot_df = prot_df.rename(columns={"Protein": "gene_name"})

    risk = "HR[95%CI]"
    if risk not in prot_df.columns: risk = "OR[95%CI]"
    risk_sm = risk.split("[")[0]
    
    #check if x is string before splitting
    prot_df[risk_sm] = prot_df[risk].apply(lambda x: float(x.split(" ")[0]) if isinstance(x, str) else x)
    prot_df[f"log{risk_sm}"] = np.log(prot_df[risk_sm])

    prot_spec_id = pd.merge(left=prot_df, right=gene_names, on="gene_name", how="left")
    
    a = prot_spec_id[prot_spec_id.duplicated('gene', keep=False)].sort_values(by="gene_name")
    pair = [item for sublist in [i.split("_") for i in a["gene_name"].tolist() if "_" in i] for item in sublist]
    
    prot_spec_id['gene'] = prot_spec_id.apply(lambda row: mappings.get(row['gene_name'], row['gene']), axis=1)
    
    gene_name_pairs = gene_names[gene_names["gene_name"].isin(pair)].drop_duplicates(subset="gene_name").rename(columns={"gene_name": "gene_name_single"})
    prot_spec_id["gene_name_single"] = prot_spec_id["gene_name"].apply(lambda x: x.split("_")[0] if "_" in x else x)
    prot_spec_id = prot_spec_id.merge(gene_name_pairs, on='gene_name_single', how='left', suffixes=('', '_small'))
    prot_spec_id['gene'] = prot_spec_id['gene'].fillna(prot_spec_id['gene_small'])
    prot_spec_id.drop(columns=['gene_small', "gene_name_single"], inplace=True)

    larger_pval = prot_spec_id[prot_spec_id["gene_name"].isin(["NPPB", "NTproBNP"])]["P_value"].max()
    dropped_prot_rows = prot_spec_id[(prot_spec_id["gene_name"].isin(["NPPB", "NTproBNP"])) & (prot_spec_id["P_value"] == larger_pval)]
    if not dropped_prot_rows.empty:
        dropped_prot = dropped_prot_rows["gene_name"].item()
        prot_spec_id = prot_spec_id[prot_spec_id["gene_name"] != dropped_prot]
    
    prot_spec_id = prot_spec_id.drop_duplicates(subset="gene")
    prot_uniq_genes = set(prot_spec_id["gene"].tolist()) - set(atlas.index.tolist())
    return prot_spec_id[~prot_spec_id["gene"].isin(list(prot_uniq_genes))].drop_duplicates(subset="gene", keep="first")

def remove_na_from_training_data(X_df: pd.DataFrame, prot_spec_final: pd.DataFrame):
    na_genes = X_df[X_df.isna().any(axis=1)].index.tolist()
    if len(na_genes) > 0:
        X_df = X_df[~X_df.index.isin(na_genes)].copy()
        prot_spec_final = prot_spec_final[~prot_spec_final["gene"].isin(na_genes)].copy()
    return X_df, prot_spec_final

def prep_data(args, df: pd.DataFrame, col: str):
    df["-log10(pval)"] = -np.log10(df["P_value"])
    df["z_score"] = (2*(df[col] > 1) - 1) * df["P_value"].apply(lambda x: stats.norm.isf(x / 2))
    max_non_inf = df.loc[df["-log10(pval)"] != np.inf, "-log10(pval)"].max()
    df = df.replace([np.inf, -np.inf], max_non_inf)
    df["-log10(pval)_minmax"] = (df["-log10(pval)"] - df["-log10(pval)"].min()) / (df["-log10(pval)"].max() - df["-log10(pval)"].min())
    y = df[args.output_label]
    if args.abs_hr: y = np.abs(y)
    weight_col = "-log10(pval)_minmax" if args.gene_weight_minmax else "-log10(pval)"
    return df, y, weight_col

# VM adaptation

LSF_FLAGS_WITH_VALUE = {"-P", "-q", "-n", "-W", "-R", "-M", "-J", "-oo", "-eo", "-o", "-e", "-cwd", "-L"}
LSF_FLAGS_NO_VALUE = {"-Is", "-I", "-K", "-B", "-N"}

def _have_lsf() -> bool:
    return shutil.which("bsub") is not None

def _strip_bsub_flags(cmd: List[str]):
    """
    Extracts core command and the paths for stdout/stderr logs.
    Works for both single-string commands and list commands.
    """
    if not cmd: return [], None, None
    
    out_log, err_log = None, None
    cleaned = []
    
    # If the command is a list, process it
    # Drop 'bsub' if it's the first element
    work_cmd = cmd[1:] if cmd[0] == "bsub" else cmd
    
    i = 0
    while i < len(work_cmd):
        tok = work_cmd[i]
        # Look for log flags and capture the next element (the path)
        if tok in ["-oo", "-o"]:
            out_log = work_cmd[i+1]
            i += 2
        elif tok in ["-eo", "-e"]:
            err_log = work_cmd[i+1]
            i += 2
        # Skip other LSF specific flags that take a value
        elif tok in ["-P", "-q", "-n", "-W", "-R", "-M", "-J", "-L", "-cwd"]:
            i += 2
        # Skip LSF specific flags that don't take a value
        elif tok in ["-Is", "-I", "-K", "-B", "-N"]:
            i += 1
        else:
            cleaned.append(tok)
            i += 1
            
    # Clean up double bash calls if they exist
    if len(cleaned) >= 2 and cleaned[0].endswith("bash") and cleaned[1] == "bash":
        cleaned = cleaned[1:]
        
    return cleaned, out_log, err_log

def submit_job_and_wait(bsub_cmd: List[str], wait_time: int = 10) -> None:
    if not bsub_cmd: return
    
    # Check if we are on HPC or VM
    have_lsf = shutil.which("bsub") is not None
    
    if not have_lsf:
        #VM MODE: Run locally and create log files
        local_cmd, out_log, err_log = _strip_bsub_flags(bsub_cmd)
        
        logging.error("VM Mode: Running locally. Outputting to: %s", out_log)

        # Execute the command and capture output to memory
        res = subprocess.run(local_cmd, capture_output=True, text=True)

        # MANUALLY write the .stdout file
        if out_log:
            os.makedirs(os.path.dirname(out_log), exist_ok=True)
            with open(out_log, "w") as f:
                f.write(res.stdout)
        
        # MANUALLY write the .stderr file
        if err_log:
            os.makedirs(os.path.dirname(err_log), exist_ok=True)
            with open(err_log, "w") as f:
                f.write(res.stderr)

        if res.returncode != 0:
           
            print(f"Error in {local_cmd[-1]}: {res.stderr}")
            raise RuntimeError(f"Local command failed (exit={res.returncode}). Logs saved.")
        return

    # HPC MODE: Use standard bsub/bjobs logic

    logging.error("HPC Mode: Submitting to LSF...")
    result = subprocess.run(bsub_cmd, capture_output=True, text=True)
    m = re.search(r"<(\d+)>", result.stdout)
    if not m:
        raise RuntimeError(f"bsub failed: {result.stderr}")
    
    jobid = m.group(1)
    time.sleep(5)
    while True:
        bjobs = subprocess.run(["bjobs", "-noheader", "-o", "stat", jobid], capture_output=True, text=True)
        stat = bjobs.stdout.strip()
        if stat in ("DONE", "EXIT"): break
        time.sleep(wait_time)