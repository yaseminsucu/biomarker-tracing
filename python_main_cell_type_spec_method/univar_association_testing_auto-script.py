import os, subprocess, argparse
from utils import submit_job_and_wait
# Instruction to run
# python $proj_path/scripts/univar_association_testing_auto-script.py --atlas_path $proj_path/results/atlas_data/atlas_all_cell_tissues.tsv --save_path_suffix dep-var-z_mean-covar-F_z-trans --z_transform 1
# proj_path="/sc/arion/projects/DiseaseGeneCell/Huang_lab_project/BioResNetwork/Phuc/projects/Alzheimer/human_atlas/sub_projects/plasma_proteome"


BASH_SCRIPT_DIR = "/mnt/vm-shared-storage/biomarker-tracing/bash_scripts/univar_association_testing.sh"
DISEASE_PROT_DIR = "/mnt/vm-shared-storage/plasma_proteome_gnpc"


def main(in_args):
    try:
        print(f"[DEBUG] Starting main()")
        print(f"[DEBUG] disease_name: {in_args.disease_name}")
        print(f"[DEBUG] disease_prot_dir: {in_args.disease_prot_dir}")
        
        base_path = f"{in_args.disease_prot_dir}/{in_args.disease_type}_popu-{in_args.popu_type}"
        print(f"[DEBUG] base_path: {base_path}")
        print(f"[DEBUG] base_path exists: {os.path.exists(base_path)}")
        
        if not os.path.exists(base_path):
            print(f"[ERROR] base_path does not exist!")
            return
        
        diseases = os.listdir(base_path)
        print(f"[DEBUG] diseases found: {diseases}")
        
        if in_args.save_path_suffix != "":
            in_args.save_path_suffix = f"_{in_args.save_path_suffix}"
        
        save_path = f"{in_args.save_path}/{in_args.disease_type}_popu-{in_args.popu_type}/univar_association_testing{in_args.save_path_suffix}"
        print(f"[DEBUG] save_path: {save_path}")
        
        lsf_params = ["-J", "atlas", "-P", "acc_DiseaseGeneCell", "-n", "1", "-W", "1:30", "-R", "rusage[mem=2000]", "-M", "20000", "-L", "/bin/bash"]

        for disease in sorted(diseases):
            print(f"[DEBUG] Processing disease file: {disease}")
            
            if len(disease.split(".")) == 2:
                dis_name = disease.split(".")[0]
            else:
                dis_name = ".".join(disease.split(".")[:-1])

            dis_name = disease.split(".")[0]
            print(f"[DEBUG] dis_name: {dis_name}, comparing to: {in_args.disease_name}")
            
            save_full_path = os.path.join(save_path, dis_name)
            os.makedirs(save_full_path, exist_ok=True)
            log_path = save_full_path

            args = [
                in_args.atlas_path, in_args.atlas_smal_path, base_path,
                save_full_path, dis_name, str(in_args.abs_hr), in_args.output_label,
                in_args.covar_df, str(in_args.covar_gini), str(in_args.ztransform_type)
            ]
            command = ["bsub"] + lsf_params + ["-oo", f"{log_path}/{disease}.stdout", "-eo", f"{log_path}/{disease}.stderr"] + ["bash", BASH_SCRIPT_DIR] + args
            
            print(f"[DEBUG] Checking condition: dis_name({dis_name}) == in_args.disease_name({in_args.disease_name})")
            
            if dis_name == in_args.disease_name:
                print(f"[DEBUG] Condition matched! Submitting bsub job and waiting...")
                submit_job_and_wait(command, wait_time=10)
                print(f"[DEBUG] Job completed!")
                break
            elif in_args.disease_name == "":
                print(f"[DEBUG] disease_name is empty, submitting all")
                submit_job_and_wait(command, wait_time=10)                
    except Exception as e:
        print(f"[ERROR] Exception in main: {e}")
        import traceback
        traceback.print_exc()



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("--atlas_smal_path", required=True, type=str)
    parser.add_argument("--atlas_path", required=True, type=str)
    parser.add_argument("--save_path", required=True, type=str)
    parser.add_argument("--save_path_suffix", required=True, type=str, default="")

    parser.add_argument("--disease_prot_dir", required=False, default=DISEASE_PROT_DIR)
    parser.add_argument("--bash_script_dir", required=False, default=BASH_SCRIPT_DIR)
    parser.add_argument("--output_label", required=False, type=str, default="HR")

    parser.add_argument("--disease_name", required=False, type=str, default="")    
    parser.add_argument("--disease_type", required=False, type=str, default="incident")
    parser.add_argument("--popu_type", required=False, type=str, default="all")

    parser.add_argument("--covar_df", type=str, default="None", help="Covariate df")
    parser.add_argument("--covar_gini", type=int, default=0, help="Add Gini coefficient as covariate")
    parser.add_argument("--abs_hr", type=int, default=0, help="Whether to set HR to abs value")
    parser.add_argument("--ztransform_type", type=int, default=1, help="Whether to z transform on each cell type (1) or each gene (2)")
    
    args = parser.parse_args()
    main(args)

