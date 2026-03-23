#!/usr/bin/env python3
"""Generate the GitHub social preview PNG (1280x640) from kagen_logo.svg."""
import os, subprocess
from PIL import Image

DIR = os.path.dirname(os.path.abspath(__file__))
LOGO_SVG = os.path.join(DIR, "..", "kagen_logo.svg")
TEMP_PNG = os.path.join(DIR, "logo-temp.png")
OUT_PNG = os.path.join(DIR, "social-preview.png")

subprocess.run(["inkscape", LOGO_SVG, "--export-type=png",
    f"--export-filename={TEMP_PNG}", "--export-height=600",
    "--export-background-opacity=0"], check=True, capture_output=True)

canvas = Image.new("RGBA", (1280, 640), (255, 255, 255, 255))
logo = Image.open(TEMP_PNG).convert("RGBA")
x = (1280 - logo.width) // 2
y = (640 - logo.height) // 2
canvas.paste(logo, (x, y), logo)
canvas.convert("RGB").save(OUT_PNG, quality=95)
os.remove(TEMP_PNG)
print(f"Saved {OUT_PNG}")
