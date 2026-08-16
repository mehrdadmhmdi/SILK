#!/usr/bin/env python3
"""
Publication-quality MACS figures, rebuilt from the result CSVs.

Replaces the cramped ggplot output (truncated legend, oversized fonts,
overlapping tick labels). Produces readable, distinguishable series with a
clean shared legend and no redundant in-figure titles (the LaTeX captions
carry the description).

Usage:
    python3 make_macs_figures.py
Reads:  results/macs_metrics_per_horizon.csv and macs_predictions.csv
Writes: figures/macs_*  (.png and .pdf), matching the manuscript figure keys
"""
import os
try:
    import numpy as np
    import pandas as pd
    import matplotlib
except ModuleNotFoundError as exc:
    raise SystemExit(
        "The optional Python renderer requires numpy, pandas, and matplotlib. "
        "The authoritative R pipeline (02_results.R) does not require it."
    ) from exc
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get(
    "MACS_RESULTS_DIR", os.environ.get("MACS_RES", os.path.join(HERE, "results"))
)
FIG = os.environ.get(
    "MACS_FIGURES_DIR", os.environ.get("MACS_FIGOUT", os.path.join(HERE, "figures"))
)
os.makedirs(FIG, exist_ok=True)
DPI = 400

# ---- house style ------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": DPI,
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 9.5,
    "axes.labelsize": 9,
    "axes.edgecolor": "#555555",
    "axes.linewidth": 0.7,
    "axes.grid": True,
    "grid.color": "#E4E4E4",
    "grid.linewidth": 0.5,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 8.5,
    "legend.handlelength": 2.2,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
})

UIUC_ORANGE = "#FF5F05"
UIUC_BLUE = "#13294B"

# display label -> (color, linestyle, marker)
STYLE = {
    "SILK-Cox (Gaussian)":                (UIUC_ORANGE, "-",  "o"),
    "Parametric MMLM-Cox":                (UIUC_BLUE,   "-.", "^"),
    "Recorded-clock Cox":                 ("#7F7F7F",   (0, (1, 1)), "x"),
    "Recorded-clock Cox (same features)": ("#B0B0B0",   (0, (3, 1)), "s"),
}
ORDER = list(STYLE.keys())
LSMAP = {  # raw method id -> display label
    "SILK-Cox": "SILK-Cox (Gaussian)",
    "MMLM-Recorded": "Parametric MMLM-Cox",
    # Primary benchmark: dynamic-prediction Cox using only the recorded time
    # information (with origin error) plus baseline demographics.
    "Cox-Recorded": "Recorded-clock Cox",
    # Optional single-channel ablation sharing SILK's full predictor map.
    "Cox-Recorded-SameFeature": "Recorded-clock Cox (same features)",
}
FACET_TITLE = {"fixed_1y": "1-year landmark", "fixed_2y": "2-year landmark",
               "fixed_3y": "3-year landmark", "fixed_5y": "5-year landmark"}


def _legend_handles(labels=None):
    labels = ORDER if labels is None else labels
    return [Line2D([0], [0], color=STYLE[l][0], linestyle=STYLE[l][1],
                   marker=STYLE[l][2], markersize=5.5, linewidth=1.6,
                   markeredgecolor=STYLE[l][0],
                   markerfacecolor="none" if STYLE[l][2] in ("x",) else STYLE[l][0])
            for l in labels]


def _load_metrics():
    m = pd.read_csv(os.path.join(RES, "macs_metrics_per_horizon.csv"))
    m["label"] = m["method"].map(LSMAP)
    return m


def _facet_grid(metric, se_col, ylab, outfile, ylim=None):
    m = _load_metrics()
    sets = ["fixed_1y", "fixed_2y", "fixed_3y", "fixed_5y"]
    horizons = sorted(m["horizon"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.2), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, ls in zip(axes, sets):
        sub = m[m["landmark_set"] == ls]
        for lab in ORDER:
            d = sub[sub["label"] == lab].sort_values("horizon")
            if d.empty:
                continue
            color, style, marker = STYLE[lab]
            mfc = "none" if marker == "x" else color
            ax.plot(d["horizon"], d[metric], linestyle=style, color=color,
                    marker=marker, markersize=5, linewidth=1.5,
                    markerfacecolor=mfc, markeredgecolor=color,
                    label=lab, zorder=3, clip_on=True)
            if se_col and se_col in d:
                mean_col = metric.replace("pooled", "mean")
                ax.errorbar(d["horizon"], d[mean_col], yerr=1.96 * d[se_col],
                            fmt="none", ecolor=color, elinewidth=0.8,
                            capsize=2, alpha=0.5, zorder=2)
        ax.set_title(FACET_TITLE[ls], fontweight="bold", fontsize=9)
        ax.set_xticks(horizons)
        if ylim:
            ax.set_ylim(*ylim)
        ax.margins(x=0.08)
    for ax in axes[2:]:
        ax.set_xlabel("Prediction horizon (years)")
    for ax in (axes[0], axes[2]):
        ax.set_ylabel(ylab)
    present = [lab for lab in ORDER if (m["label"] == lab).any()]
    fig.legend(_legend_handles(present), present, loc="lower center",
               ncol=min(3, max(1, len(present))), frameon=False,
               bbox_to_anchor=(0.5, -0.06), columnspacing=1.4)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(os.path.join(FIG, outfile + ".png"))
    fig.savefig(os.path.join(FIG, outfile + ".pdf"))
    plt.close(fig)
    print("wrote", outfile)


def fig1_auc():
    # Y-limits derive from the data so no pooled AUC can fall off the plot.
    _facet_grid("auc_pooled", None, "Time-dependent IPCW AUC",
                "macs_auc_by_horizon")


def fig2_brier():
    _facet_grid("brier_pooled", None, "IPCW Brier score",
                "macs_brier_by_horizon")


def _km_censor(u, d):
    """Return landmark-specific censoring survival G(t) and G(t-)."""
    uniq = np.unique(u)
    km = 1.0
    keys, vals = [], []
    for t in uniq:
        dt = np.sum((u == t) & (d == 0))   # censoring "events"
        nt = np.sum(u >= t)
        if nt > 0:
            km *= (1 - dt / nt)
        keys.append(t); vals.append(km)
    keys = np.array(keys); vals = np.array(vals)

    def G(x, left_limit=False):
        x = np.atleast_1d(np.asarray(x, float))
        out = np.ones(len(x))
        for i, xi in enumerate(x):
            side = "left" if left_limit else "right"
            index = np.searchsorted(keys, xi, side=side) - 1
            out[i] = vals[index] if index >= 0 else 1.0
        return np.clip(out, 1e-6, 1.0)
    return G


def fig3_calibration():
    """Landmark-specific IPCW calibration at the prespecified 5-year horizon."""
    pred = pd.read_csv(os.path.join(RES, "macs_predictions.csv"))
    h = 5 if 5 in pred["horizon"].unique() else pred["horizon"].max()
    sets = [x for x in ["fixed_1y", "fixed_2y", "fixed_3y", "fixed_5y"]
            if x in pred["landmark_set"].unique()]
    if not sets:
        print("skipped fig3_calibration: no fixed-landmark predictions")
        return
    fig, axes = plt.subplots(2, 2, figsize=(7.0, 5.7), sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, landmark_set in zip(axes, sets):
        landmark = pred[pred["landmark_set"] == landmark_set]
        outc = landmark.drop_duplicates("subject_id")
        G = _km_censor(outc["residual_time"].values.astype(float),
                       outc["event_status"].values.astype(int))
        for lab in ORDER:
            raw = [k for k, v in LSMAP.items() if v == lab][0]
            z = landmark[(landmark["method"] == raw) &
                         (landmark["horizon"] == h)].copy()
            if len(z) < 40:
                continue
            U = z["residual_time"].values.astype(float)
            D = z["event_status"].values.astype(int)
            y = np.where((D == 1) & (U <= h), 1.0, np.where(U > h, 0.0, np.nan))
            w = np.where((D == 1) & (U <= h), 1.0 / G(U, left_limit=True),
                         np.where(U > h, 1.0 / G(np.full(len(U), h)), 0.0))
            p = z["risk_pred"].values.astype(float)
            keep = np.isfinite(y) & (w > 0)
            y, w, p = y[keep], w[keep], p[keep]
            if len(y) < 40:
                continue
            q = np.unique(np.quantile(p, np.linspace(0, 1, 6)))
            if len(q) < 3:
                continue
            b = np.clip(np.digitize(p, q[1:-1]), 0, len(q) - 2)
            xs, ys = [], []
            for bi in np.unique(b):
                mask = b == bi
                xs.append(np.average(p[mask], weights=w[mask]))
                ys.append(np.average(y[mask], weights=w[mask]))
            color, style, marker = STYLE[lab]
            mfc = "none" if marker == "x" else color
            ax.plot(xs, ys, linestyle=style, color=color, marker=marker,
                    markersize=5, linewidth=1.4, markerfacecolor=mfc,
                    markeredgecolor=color, label=lab, zorder=3)
        ax.plot([0, 1], [0, 1], "--", color="#999999", linewidth=1, zorder=1)
        ax.set_aspect("equal")
        ax.set_title(FACET_TITLE.get(landmark_set, landmark_set),
                     fontweight="bold", fontsize=9)
    for ax in axes[len(sets):]:
        ax.set_visible(False)
    visible_axes = axes[:len(sets)]
    finite_max = 0.05
    for ax in visible_axes:
        for line in ax.lines[:-1]:
            if len(line.get_xdata()):
                finite_max = max(finite_max, np.nanmax(line.get_xdata()),
                                 np.nanmax(line.get_ydata()))
    axis_max = min(1.0, max(0.10, 1.05 * finite_max))
    for index, ax in enumerate(visible_axes):
        ax.set_xlim(0, axis_max)
        ax.set_ylim(0, axis_max)
        if index >= 2:
            ax.set_xlabel("Predicted risk")
        if index % 2 == 0:
            ax.set_ylabel("Observed risk (IPCW)")
    present = [lab for lab in ORDER
               if any(lab == line.get_label() for ax in visible_axes
                      for line in ax.lines)]
    if not present:
        present = ORDER[:3]
    fig.legend(_legend_handles(present), present, loc="lower center", ncol=2,
               frameon=False, bbox_to_anchor=(0.5, -0.02), columnspacing=1.4)
    fig.suptitle(f"{h:g}-year risk calibration", fontsize=10, fontweight="bold")
    fig.tight_layout(rect=[0, 0.08, 1, 1])
    fig.savefig(os.path.join(FIG, "macs_calibration.png"))
    fig.savefig(os.path.join(FIG, "macs_calibration.pdf"))
    plt.close(fig)
    print("wrote macs_calibration")


def fig4_km():
    """KM curves by SILK (Gaussian) 5-year risk tertile at the 5-year landmark.

    The 5-year landmark matches the R renderer (02_results.R) so both
    pipelines draw the same figure under the same filename.
    """
    pred = pd.read_csv(os.path.join(RES, "macs_predictions.csv"))
    z = pred[(pred["method"] == "SILK-Cox") & (pred["horizon"] == 5) &
             (pred["landmark_set"] == "fixed_5y")].copy()
    if len(z) < 30:
        z = pred[(pred["method"] == "SILK-Cox") & (pred["horizon"] == 5) &
                 (pred["landmark_set"] == "fixed_3y")].copy()
    if len(z) < 30:
        z = pred[(pred["method"] == "SILK-Cox") & (pred["horizon"] == 5)].copy()
    if len(z) < 30:
        print("skipped fig4_km_risk_strata: fewer than 30 SILK-Cox predictions")
        return
    z = z.drop_duplicates("subject_id")
    U = z["residual_time"].values.astype(float)
    D = z["event_status"].values.astype(int)
    p = z["risk_pred"].values.astype(float)
    ter = np.quantile(p, [1 / 3, 2 / 3])
    grp = np.where(p <= ter[0], 0, np.where(p <= ter[1], 1, 2))
    names = ["Low risk (bottom third)", "Medium risk", "High risk (top third)"]
    styles = [(UIUC_BLUE, ":", 1.5),
              (UIUC_BLUE, "--", 1.8),
              (UIUC_BLUE, "-", 2.2)]

    def km(u, d):
        uniq = np.unique(u[d == 1])
        t = [0.0]; s = [1.0]; surv = 1.0
        for tt in uniq:
            nt = np.sum(u >= tt); dt = np.sum((u == tt) & (d == 1))
            if nt > 0:
                surv *= (1 - dt / nt)
            t.append(tt); s.append(surv)
        return np.array(t), np.array(s)

    fig, ax = plt.subplots(figsize=(5.6, 4.2))
    xmax = min(10, float(np.nanmax(U)))
    for g, nm, (color, linestyle, linewidth) in zip([0, 1, 2], names, styles):
        mask = grp == g
        t, s = km(U[mask], D[mask])
        ax.step(t, s, where="post", color=color, linestyle=linestyle,
                linewidth=linewidth, label=nm)
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time from landmark (years)")
    ax.set_ylabel("AIDS-free survival probability")
    ax.legend(loc="lower left", frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "macs_km_risk_strata.png"))
    fig.savefig(os.path.join(FIG, "macs_km_risk_strata.pdf"))
    plt.close(fig)
    print("wrote fig4_km_risk_strata")


if __name__ == "__main__":
    fig1_auc()
    fig2_brier()
    fig3_calibration()
    fig4_km()
    print("All MACS figures regenerated at %d dpi -> %s" % (DPI, FIG))
