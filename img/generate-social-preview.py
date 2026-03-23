#!/usr/bin/env python3
"""Generate the GitHub social preview PNG (1280x640) from social-preview.svg."""
import os, subprocess
DIR = os.path.dirname(os.path.abspath(__file__))
subprocess.run(["inkscape", os.path.join(DIR, "social-preview.svg"), "--export-type=png", f"--export-filename={os.path.join(DIR, 'social-preview.png')}", "--export-width=1280", "--export-height=640"], check=True, capture_output=True)
