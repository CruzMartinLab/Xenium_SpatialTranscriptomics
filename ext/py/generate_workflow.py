#!/usr/bin/env python3
"""
generate_sbatch.py: create both a parameterized Makefile and a single-SLURM-job submission script.

Usage:
  python generate_sbatch.py \
      --config config.json \
      --output submit.sh

Requires a Makefile template (default: Makefile.template) with Python-style placeholders,
like:

  DATADIR = {DATADIR}
  THREADS = {THREADS}
  HEXGRIDS = {HEXGRIDS}
  TOPICS  = {TOPICS}
  # etc.

This script will:
 1. Load JSON config with "job" and "workflow".
 2. Fill in the Makefile template and write to Makefile.generated.
 3. Emit an sbatch wrapper (submit.sh) that calls
      make -f Makefile.generated -j {THREADS}

Once generated, you can:
  sbatch submit.sh
"""
import argparse
import json
import os
import sys


def render_makefile(template_path, out_path, wf):
    # Prepare formatting dict: uppercase keys, lists to space-joined strings
    fmt = {}
    for k, v in wf.items():
        key = k.upper()
        if isinstance(v, list):
            fmt[key] = " ".join(str(x) for x in v)
        else:
            fmt[key] = str(v)
    # Read template
    try:
        with open(template_path) as tf:
            tpl = tf.read()
    except IOError:
        sys.exit(f"Error: cannot read Makefile template '{template_path}'")
    # Render
    try:
        rendered = tpl.format(**fmt)
    except KeyError as e:
        sys.exit(f"Error: placeholder {e} not found in workflow config.")
    # Write out
    with open(out_path, 'w') as of:
        of.write(rendered)
    print(f"Generated Makefile: {out_path}")


def write_sbatch_script(job, makefile_out, sbatch_out, pyenv=None):
    lines = ['#!/bin/bash']
    # SBATCH headers
    for key, val in job.items():
        if key != "extra_lines":
            lines.append(f"#SBATCH --{key}={val}")
        else:
            lines.append(val)
    lines.append("")
    if pyenv:
        lines.append(f"# Activate Python environment: {pyenv}")
        lines.append(f"source {pyenv}/bin/activate")
    lines.append("\nset -euo pipefail")
    lines.append("")
    lines.append(f"# Run the generated Makefile")
    lines.append(f"make -f {makefile_out}")
    lines.append("")
    with open(sbatch_out, 'w') as f:
        f.write("\n".join(lines))
    os.chmod(sbatch_out, 0o755)
    print(f"Generated SBATCH wrapper: {sbatch_out}")


def main():
    p = argparse.ArgumentParser(description="Generate Makefile + SLURM script from config")
    p.add_argument('-c', '--config', required=True, help="JSON config file path")
    p.add_argument('-o', '--output', default="", help="Output sbatch script path")
    p.add_argument('-t', '--template', required=True,
                   help="Makefile template path (with {PLACEHOLDERS})")
    p.add_argument('-m', '--makefile', default='Makefile',
                   help="Output Makefile path")
    args = p.parse_args()

    # Load config
    cfg = json.load(open(args.config))
    job = cfg.get('job', {})
    wf  = cfg.get('workflow', {})
    pyenv = cfg.get('pyenv', None)
    # Render Makefile
    render_makefile(args.template, args.makefile, wf)
    # Emit sbatch wrapper
    if (args.output):
        write_sbatch_script(job, args.makefile, args.output, pyenv)
    # check temporary directory
    tmpdir= wf.get("tmpdir")
    if os.path.exists(tmpdir):
        print(f"Please removed temporary directory or all its content before running the workflow: {tmpdir}")

if __name__ == '__main__':
    main()
