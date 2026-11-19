"""Shared plotting helpers: Times New Roman (with serif fallback), 600-dpi JPEG
export, data/metric loaders, deterministic layout, and colour maps."""
import csv, json, os, io, glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from PIL import Image

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA = os.path.join(HERE, "data"); OUT = os.path.join(HERE, "outputs")
FIGS = os.path.join(HERE, "figures")

# Times New Roman if present (e.g. ~/.fonts or figures/fonts), else default serif.
for pat in ("/root/.fonts/times*.ttf", os.path.join(FIGS, "fonts", "*.ttf")):
    for f in glob.glob(pat):
        try: fm.fontManager.addfont(f)
        except Exception: pass
_FAM = "Times New Roman" if any("Times New Roman" == f.name for f in fm.fontManager.ttflist) else "serif"
plt.rcParams.update({"font.family": _FAM, "font.size": 11})

CCOL = {"1": "#2e6f9e", "2": "#2e8b57", "3": "#c0392b", "4": "#8e44ad", "5": "#e67e22"}
ECOL = {"cross": "#b03a2e", "intra_predator": "#2e8b57", "intra_prey": "#1f4e79"}
BB_COL, ST_COL = "#2e8b57", "#1f4e79"


def load_nodes(): return list(csv.DictReader(open(os.path.join(DATA, "ppi_nodes.csv"))))
def load_edges(): return list(csv.DictReader(open(os.path.join(DATA, "ppi_edges.csv"))))
def load_metrics(): return json.load(open(os.path.join(OUT, "metrics.json")))
def load_layout():
    p = os.path.join(DATA, "network_layout.json")
    return {k: tuple(v) for k, v in json.load(open(p)).items()} if os.path.exists(p) else None


def savejpg(fig, name, dpi=600):
    os.makedirs(FIGS, exist_ok=True)
    buf = io.BytesIO()
    fig.savefig(buf, dpi=dpi, format="png", facecolor="white",
                bbox_inches="tight", pad_inches=0.15)
    buf.seek(0)
    Image.open(buf).convert("RGB").save(os.path.join(FIGS, f"{name}.jpg"),
                                        dpi=(dpi, dpi), quality=95)
    plt.close(fig)
    print("saved", name + ".jpg")
