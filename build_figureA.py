"""
build_figureA.py — rebuild the "Woody Biomass Figure" (Claude Design export) as a clean, self-
contained static HTML figure, resolving the Claude Design templating:
  - the JS-computed annual tree-ring radial gradient (ported from the export's DCLogic.buildRings)
  - <sc-if> conditionals -> their default values (rings on, linkage notes on, pie not donut)
  - <x-dc>/support.js wrapper removed

Writes docs/fig/figureA.html (+ copies the high-res molecule renders into docs/fig/assets/).
That static file is rasterized to docs/img/figureA-woody-biomass.png and embedded in learn.html.

    python build_figureA.py
"""
import os
import shutil

ROOT = os.path.dirname(os.path.abspath(__file__))
FIG_DIR = os.path.join(ROOT, "docs", "fig")
ASSET_DIR = os.path.join(FIG_DIR, "assets")
EXPORT_DIR = os.path.join(ROOT, "export")

IMAGES = ["lignin_trimer", "lignin_b5", "coniferyl", "ggm_extended",
          "agx_compact", "cellohexaose", "cellooctaose"]


def build_ring_gradient():
    """Port of the export's DCLogic.buildRings(): concentric earlywood/latewood annual rings."""
    R, r, w, f, minW = 171, 7, 22, 0.88, 4
    early = "rgba(170,120,60,0)"
    earlyHi = "rgba(255,243,216,0.16)"
    late = "rgba(52,29,11,0.58)"
    pct = lambda v: f"{v / R * 100:.2f}%"
    stops = ["rgba(58,36,16,0.9) 0%", "rgba(58,36,16,0.9) " + pct(r)]
    while r < R:
        end = min(R, r + w)
        span = end - r
        stops.append(early + " " + pct(r))
        stops.append(earlyHi + " " + pct(r + span * 0.28))
        stops.append(early + " " + pct(r + span * 0.6))
        stops.append(late + " " + pct(end))
        r = end
        w = max(minW, w * f)
    return "radial-gradient(circle at 50% 50%, " + ", ".join(stops) + ")"


FIGURE = """<!DOCTYPE html>
<html><head><meta charset="utf-8">
<style>*{{box-sizing:border-box;margin:0;padding:0;}}
body{{background:#ffffff;font-family:'Helvetica Neue',Helvetica,Arial,sans-serif;}}</style></head>
<body>
<div id="figcard" style="width:1660px;background:#ffffff;padding:30px 30px 34px;">
  <div style="font-size:30px;line-height:1.4;color:#1a1a1a;margin:0 4px 22px;width:1540px;">
    <span style="font-weight:700;">Figure A.</span> Molecular diagrams of cellulose, hemicellulose, and lignin&mdash;the most abundant molecules in woody biomass.
  </div>
  <div style="position:relative;width:1600px;height:1180px;border:1px solid #c9c6bf;background:#ffffff;overflow:hidden;">
    <svg viewBox="0 0 1600 1180" width="1600" height="1180" style="position:absolute;top:0;left:0;pointer-events:none;z-index:1;">
      <defs><marker id="ah" markerWidth="10" markerHeight="10" refX="6.5" refY="3.4" orient="auto">
        <path d="M0,0 L7.4,3.4 L0,6.8 Z" fill="#8a6a4a"></path></marker></defs>
      <path d="M684,422 L516,418" stroke="#8a6a4a" stroke-width="2" fill="none" marker-end="url(#ah)"></path>
      <path d="M946,454 L1086,418" stroke="#8a6a4a" stroke-width="2" fill="none" marker-end="url(#ah)"></path>
      <path d="M756,736 L800,832" stroke="#8a6a4a" stroke-width="2" fill="none" marker-end="url(#ah)"></path>
    </svg>
    <div style="position:absolute;left:598px;top:296px;width:404px;text-align:center;z-index:2;">
      <div style="font-size:22px;letter-spacing:.04em;text-transform:uppercase;color:#7a6a55;font-weight:600;">Composition of dry wood (% by mass)</div>
    </div>
    <div style="position:absolute;left:620px;top:380px;width:360px;height:360px;z-index:2;">
      <div style="position:absolute;inset:0;border-radius:50%;background:conic-gradient(from 0deg,#A9712F 0deg 107deg,#ffffff 107deg 109deg,#C68A4A 109deg 279.8deg,#ffffff 279.8deg 281.8deg,#5E3A1E 281.8deg 359deg,#ffffff 359deg 360deg);box-shadow:inset 0 0 0 4px #6b4a28, inset 0 0 0 9px #4a2f17, 0 8px 22px rgba(0,0,0,.18);"></div>
      <div style="position:absolute;inset:9px;border-radius:50%;pointer-events:none;background:{ring_gradient};"></div>
      <div style="position:absolute;inset:9px;border-radius:50%;pointer-events:none;background:radial-gradient(circle at 40% 35%, rgba(255,238,210,0.22), rgba(255,238,210,0) 46%), radial-gradient(circle at 62% 66%, rgba(34,17,4,0.30), rgba(34,17,4,0) 60%);"></div>
      <div style="position:absolute;left:50%;top:50%;width:13px;height:13px;border-radius:50%;background:#3c2611;transform:translate(-50%,-50%);box-shadow:0 0 6px rgba(0,0,0,.35);"></div>
      <div style="position:absolute;left:258px;top:118px;transform:translate(-50%,-50%);text-align:center;color:#fff;text-shadow:0 1px 4px rgba(0,0,0,.55);">
        <div style="font-size:26px;font-weight:700;">Hemicellulose</div><div style="font-size:21px;font-weight:500;">23&ndash;32%</div></div>
      <div style="position:absolute;left:157px;top:280px;transform:translate(-50%,-50%);text-align:center;color:#2e1c0c;text-shadow:0 1px 3px rgba(255,245,225,.5);">
        <div style="font-size:29px;font-weight:700;">Cellulose</div><div style="font-size:21px;font-weight:500;">38&ndash;50%</div></div>
      <div style="position:absolute;left:116px;top:102px;transform:translate(-50%,-50%);text-align:center;color:#fff;text-shadow:0 1px 4px rgba(0,0,0,.6);">
        <div style="font-size:26px;font-weight:700;">Lignin</div><div style="font-size:21px;font-weight:500;">15&ndash;25%</div></div>
    </div>
    <div style="position:absolute;left:30px;top:110px;width:472px;height:440px;border:1px solid #ddd9d1;border-top:4px solid #5E3A1E;background:#fcfbf9;padding:16px 18px;z-index:2;">
      <div style="display:flex;align-items:baseline;gap:11px;"><span style="width:13px;height:13px;border-radius:50%;background:#5E3A1E;display:inline-block;"></span>
        <span style="font-size:33px;font-weight:700;color:#2a1c10;">Lignin</span>
        <span style="font-size:21px;font-weight:600;color:#5E3A1E;background:#f0e7dd;padding:3px 10px;border-radius:11px;">15&ndash;25%</span></div>
      <div style="font-size:21px;color:#7a6a55;margin:6px 0 8px 24px;">Amorphous phenylpropanoid network</div>
      <div style="width:100%;height:132px;"><img src="assets/lignin_trimer.png" style="width:100%;height:100%;object-fit:contain;display:block;"></div>
      <div style="font-family:Menlo,Consolas,monospace;font-size:19px;color:#8a7a64;text-align:center;">lignin substructure (&beta;-O-4 / &beta;-5)</div>
      <div style="display:flex;gap:14px;margin-top:8px;">
        <div style="flex:1;"><div style="height:112px;"><img src="assets/lignin_b5.png" style="width:100%;height:100%;object-fit:contain;display:block;"></div>
          <div style="font-family:Menlo,Consolas,monospace;font-size:19px;color:#8a7a64;text-align:center;">&beta;-5 phenylcoumaran</div></div>
        <div style="flex:1;"><div style="height:112px;"><img src="assets/coniferyl.png" style="width:100%;height:100%;object-fit:contain;display:block;"></div>
          <div style="font-family:Menlo,Consolas,monospace;font-size:19px;color:#8a7a64;text-align:center;">coniferyl alcohol &mdash; monomer</div></div></div>
      <div style="font-size:21px;font-style:italic;color:#6a5a45;margin-top:9px;">Cross-links via &beta;-O-4, &beta;-5 and &beta;&ndash;&beta; bonds.</div>
    </div>
    <div style="position:absolute;left:1098px;top:110px;width:472px;height:440px;border:1px solid #ddd9d1;border-top:4px solid #A9712F;background:#fcfbf9;padding:16px 18px;z-index:2;">
      <div style="display:flex;align-items:baseline;gap:11px;"><span style="width:13px;height:13px;border-radius:50%;background:#A9712F;display:inline-block;"></span>
        <span style="font-size:33px;font-weight:700;color:#2a1c10;">Hemicellulose</span>
        <span style="font-size:21px;font-weight:600;color:#8a5a24;background:#f4ebe0;padding:3px 10px;border-radius:11px;">23&ndash;32%</span></div>
      <div style="font-size:21px;color:#7a6a55;margin:6px 0 12px 24px;">Branched heteropolysaccharides</div>
      <div style="height:162px;"><img src="assets/ggm_extended.png" style="width:100%;height:100%;object-fit:contain;display:block;"></div>
      <div style="font-family:Menlo,Consolas,monospace;font-size:19px;color:#8a7a64;text-align:center;margin-bottom:8px;">galactoglucomannan</div>
      <div style="height:104px;"><img src="assets/agx_compact.png" style="width:100%;height:100%;object-fit:contain;display:block;"></div>
      <div style="font-family:Menlo,Consolas,monospace;font-size:19px;color:#8a7a64;text-align:center;">arabinoglucuronoxylan</div>
      <div style="font-size:21px;font-style:italic;color:#6a5a45;margin-top:9px;">&beta;-(1&rarr;4) backbone with galactose, arabinose &amp; glucuronic-acid side groups.</div>
    </div>
    <div style="position:absolute;left:360px;top:840px;width:880px;height:320px;border:1px solid #ddd9d1;border-top:4px solid #C68A4A;background:#fcfbf9;padding:16px 20px;z-index:2;">
      <div style="display:flex;align-items:baseline;gap:11px;"><span style="width:13px;height:13px;border-radius:50%;background:#C68A4A;display:inline-block;"></span>
        <span style="font-size:33px;font-weight:700;color:#2a1c10;">Cellulose</span>
        <span style="font-size:21px;font-weight:600;color:#9a6a2a;background:#f6efe4;padding:3px 10px;border-radius:11px;">38&ndash;50%</span></div>
      <div style="font-size:21px;color:#7a6a55;margin:6px 0 10px 24px;">Linear, crystalline glucan chains</div>
      <div style="display:flex;gap:28px;align-items:center;">
        <div style="flex:1;"><div style="height:160px;"><img src="assets/cellohexaose.png" style="width:100%;height:100%;object-fit:contain;display:block;"></div>
          <div style="font-family:Menlo,Consolas,monospace;font-size:19px;color:#8a7a64;text-align:center;">cellohexaose (DP 6)</div></div>
        <div style="flex:1;"><div style="height:160px;"><img src="assets/cellooctaose.png" style="width:100%;height:100%;object-fit:contain;display:block;"></div>
          <div style="font-family:Menlo,Consolas,monospace;font-size:19px;color:#8a7a64;text-align:center;">cellooctaose (DP 8)</div></div></div>
      <div style="font-size:21px;font-style:italic;color:#6a5a45;margin-top:8px;text-align:center;">Linear &beta;-(1&rarr;4)-D-glucan; hydrogen-bonds into crystalline microfibrils.</div>
    </div>
  </div>
</div>
</body></html>
"""


def main():
    os.makedirs(ASSET_DIR, exist_ok=True)
    for name in IMAGES:
        src = os.path.join(EXPORT_DIR, name + ".png")
        if not os.path.exists(src):
            raise SystemExit(f"missing high-res render {src} — run `python export_figures.py` first")
        shutil.copy2(src, os.path.join(ASSET_DIR, name + ".png"))
    html = FIGURE.format(ring_gradient=build_ring_gradient())
    with open(os.path.join(FIG_DIR, "figureA.html"), "w") as fh:
        fh.write(html)
    print(f"wrote docs/fig/figureA.html + {len(IMAGES)} assets")


if __name__ == "__main__":
    main()
