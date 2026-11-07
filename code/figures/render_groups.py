"""Render the chosen regional display and equal-scale structural panels."""
from pathlib import Path
import subprocess,sys
R=Path(__file__).resolve().parents[2]
part=sys.argv[1] if len(sys.argv)>1 else 'tu'
subprocess.run([sys.executable,str(R/'code/figures/render_regional.py'),part,'--assets',str(R/'data/array_study'),'--labels',str(R/'data/figure_source/figures'),'--windows',str(R/'data/regional_display/regional_windows_20kb_step5kb.tsv.gz'),'--out',str(R/'figures')],check=True)
