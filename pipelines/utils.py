#!/usr/bin/env python3

import os
import re
import time
import shutil
import logging
import subprocess
import argparse
import copy
import yaml
from typing import List

LSF_FLAGS_WITH_VALUE = {"-P", "-q", "-n", "-W", "-R", "-M", "-J", "-oo", "-eo", "-o", "-e", "-cwd"}
LSF_FLAGS_NO_VALUE = {"-Is", "-I", "-K", "-B", "-N"}


def _have_lsf() -> bool:
    return (shutil.which("bsub") is not None) and (shutil.which("bjobs") is not None)


def _extract_underlying_cmd_from_bsub(bsub_cmd: List[str]) -> List[str]:
    """
    Given a full bsub command list, return the underlying command list by stripping LSF flags.
    Only call this if bsub_cmd[0] == "bsub".
    """
    if not bsub_cmd:
        raise ValueError("Empty command list")
    if bsub_cmd[0] != "bsub":
        raise ValueError("Expected a bsub command")

    cmd = bsub_cmd[1:]  # drop "bsub"

    cleaned: List[str] = []
    i = 0
    while i < len(cmd):
        tok = cmd[i]
        if tok in LSF_FLAGS_WITH_VALUE:
            i += 2
            continue
        if tok in LSF_FLAGS_NO_VALUE:
            i += 1
            continue
        cleaned.append(tok)
        i += 1

    return cleaned


def submit_job_and_wait(bsub_cmd: List[str], wait_time: int = 10) -> None:
    
    if not bsub_cmd:
        raise ValueError("submit_job_and_wait: bsub_cmd is empty")

    have_lsf = _have_lsf()

    #LOCAL / VM PATH 
    if (not have_lsf) or (bsub_cmd[0] != "bsub"):
        if bsub_cmd[0] == "bsub":
            local_cmd = _extract_underlying_cmd_from_bsub(bsub_cmd)
        else:
            local_cmd = bsub_cmd  

        logging.error("Running locally:\n%s", " ".join(local_cmd))
        result = subprocess.run(local_cmd, capture_output=True, text=True)

        if result.stdout:
            logging.error("STDOUT:\n%s", result.stdout)
        if result.stderr:
            logging.error("STDERR:\n%s", result.stderr)

        if result.returncode != 0:
            raise RuntimeError(
                f"Local command failed (exit={result.returncode}).\n"
                f"CMD: {' '.join(local_cmd)}\n"
                f"STDERR (tail): {result.stderr[-2000:] if result.stderr else ''}"
            )
        return

    # HPC / LSF PATH 
    logging.error("Submitting to LSF:\n%s", " ".join(bsub_cmd))
    result = subprocess.run(bsub_cmd, capture_output=True, text=True)

    if result.returncode != 0:
        raise RuntimeError(
            f"bsub failed (exit={result.returncode}).\n"
            f"STDOUT:\n{result.stdout}\n"
            f"STDERR:\n{result.stderr}"
        )

    m = re.search(r"<(\d+)>", result.stdout)
    if not m:
        raise RuntimeError(f"Could not parse job id from bsub output:\n{result.stdout}\n{result.stderr}")

    jobid = m.group(1)
    logging.error("Submitted job %s, waiting...", jobid)

    time.sleep(10)

    while True:
        bjobs = subprocess.run(["bjobs", "-noheader", "-o", "stat", jobid],
                              capture_output=True, text=True)
        stat = bjobs.stdout.strip()

        if stat in ("DONE", "EXIT"):
            logging.error("stat=%s, job %s finished", stat, jobid)
            break

        logging.error("stat=%s, job %s still running/pending", stat, jobid)
        time.sleep(wait_time)


def create_yml_files(base_yml: str, diseases: list, save_path: str):
    os.makedirs(save_path, exist_ok=True)
    with open(base_yml, "r") as f:
        config = yaml.safe_load(f)

    for disease in sorted(diseases):
        new_config = copy.deepcopy(config)
        new_config["inputs"]["disease_name"] = [disease]

        outpath = os.path.join(save_path, f"{disease}.yml")
        with open(outpath, "w") as file:
            yaml.safe_dump(new_config, file, default_flow_style=False, sort_keys=False)


def main(args):
    diseases = [os.path.splitext(i)[0] for i in os.listdir(args.disease_path)]
    create_yml_files(args.yml_full_path, diseases, os.path.join(args.save_path, args.save_name))


if __name__ == "__main__":
    logging.basicConfig(level=logging.ERROR, format="%(levelname)s:%(message)s")

    parser = argparse.ArgumentParser()
    parser.add_argument("--save_path", type=str, required=True)
    parser.add_argument("--save_name", type=str, required=True)
    parser.add_argument("--yml_full_path", type=str, required=True)
    parser.add_argument("--disease_path", type=str, required=True)
    args = parser.parse_args()

    main(args)
