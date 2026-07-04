#!/usr/bin/env python3
"""
Publication-quality MACS figures, rebuilt from the result CSVs.

Replaces the cramped ggplot output (truncated legend, oversized fonts,
overlapping tick labels). Produces readable, distinguishable series with a
clean shared legend and no redundant in-figure titles (the LaTeX captions
carry the description).

Usage:
    python3 make_macs_figures.py
Reads:  results/macs_metrics_per_horizon.csv, macs_predictions.csv,
        macs_silk_shifts.csv
Writes: figures/fig1..fig5  (.png and .pdf)
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.environ.get("MACS_RES", os.path.join(HERE, "results"))
FIG = os.environ.get("MACS_FIGOUT", os.path.join(HERE, "figures"))
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

UIUC_ORANGE = "#E8590C"
UIUC_BLUE = "#13294B"

# display label -> (color, linestyle, marker)
STYLE = {
    "SILK (Gaussian kernel)":   (UIUC_ORANGE, "-",  "o"),
    "SILK (linear kernel)":     ("#8C564B",   "--", "s"),
    "Mixed-model landmarking":  (UIUC_BLUE,   "-.", "^"),
    "Joint model":              ("#C0392B",   ":",  "D"),
    "Landmarking":              ("#7F7F7F",   (0, (1, 1)),  "x"),
}
ORDER = list(STYLE.keys())
LSMAP = {  # raw method id -> display label
    "SILK": "SILK (Gaussian kernel)",
    "SILK-LinearMMD": "SILK (linear kernel)",
    "MMLM-Recorded": "Mixed-model landmarking",
    "JM-Recorded": "Joint model",
    "Landmark-Recorded": "Landmarking",
}
FACET_TITLE = {"fixed_1y": "1-year landmark", "fixed_2y": "2-year landmark",
               "fixed_3y": "3-year landmark", "fixed_5y": "5-year landmark"}


def _legend_handles():
    return [Line2D([0], [0], color=STYLE[l][0], linestyle=STYLE[l][1],
                   marker=STYLE[l][2], markersize=5.5, linewidth=1.6,
                   markeredgecolor=STYLE[l][0],
                   markerfacecolor="none" if STYLE[l][2] in ("x",) else STYLE[l][0])
            for l in ORDER]


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
    fig.legend(_legend_handles(), ORDER, loc="lower center", ncol=3,
               frameon=False, bbox_to_anchor=(0.5, -0.06), columnspacing=1.4)
    fig.tight_layout(rect=[0, 0.06, 1, 1])
    fig.savefig(os.path.join(FIG, outfile + ".png"))
    fig.savefig(os.path.join(FIG, outfile + ".pdf"))
    plt.close(fig)
    print("wrote", outfile)


def fig1_auc():
    _facet_grid("auc_pooled", "auc_se", "Time-dependent IPCW AUC",
                "fig1_auc_by_horizon", ylim=(0.40, 0.90))


def fig2_brier():
    _facet_grid("brier_pooled", "brier_se", "IPCW Brier score",
                "fig2_brier_by_horizon")


def _km_censor(u, d):
    """Return G(t): KM estimate of the censoring survival at times t (step, left-cont)."""
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

    def G(x):
        x = np.atleast_1d(np.asarray(x, float))
        out = np.ones(len(x))
        for i, xi in enumerate(x):
            le = keys[keys <= xi]
            out[i] = vals[np.searchsorted(keys, le[-1])] if len(le) else 1.0
        return np.clip(out, 1e-6, 1.0)
    return G


def fig3_calibration():
    """IPCW predicted-vs-observed risk (quintiles), pooled, at horizons 2 and 5 yr."""
    pred = pd.read_csv(os.path.join(RES, "macs_predictions.csv"))
    outc = pred.drop_duplicates("subject_id")
    G = _km_censor(outc["residual_time"].values.astype(float),
                   outc["event_status"].values.astype(int))
    horizons = [2, 5]
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.7))
    for ax, h in zip(axes, horizons):
        for lab in ORDER:
            raw = [k for k, v in LSMAP.items() if v == lab][0]
            z = pred[(pred["method"] == raw) & (pred["horizon"] == h)].copy()
            if len(z) < 40:
                continue
            U = z["residual_time"].values.astype(float)
            D = z["event_status"].values.astype(int)
            y = np.where((D == 1) & (U <= h), 1.0, np.where(U > h, 0.0, np.nan))
            w = np.where((D == 1) & (U <= h), 1.0 / G(U),
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
        amax = 0.6 if h == 5 else 0.35
        ax.plot([0, amax], [0, amax], "--", color="#999999", linewidth=1, zorder=1)
        ax.set_xlim(0, amax); ax.set_ylim(0, amax)
        ax.set_aspect("equal")
        ax.set_title(f"{h}-year horizon", fontweight="bold", fontsize=9)
        ax.set_xlabel("Predicted risk")
    axes[0].set_ylabel("Observed proportion (IPCW)")
    fig.legend(_legend_handles(), ORDER, loc="lower center", ncol=3,
               frameon=False, bbox_to_anchor=(0.5, -0.10), columnspacing=1.4)
    fig.tight_layout(rect=[0, 0.04, 1, 1])
    fig.savefig(os.path.join(FIG, "fig3_calibration.png"))
    fig.savefig(os.path.join(FIG, "fig3_calibration.pdf"))
    plt.close(fig)
    print("wrote fig3_calibration")


def fig4_shifts():
    s = pd.read_csv(os.path.join(RES, "macs_silk_shifts.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.5))
    ax = axes[0]
    ax.hist(s["e_hat"], bins=40, color=UIUC_ORANGE, edgecolor="white",
            linewidth=0.4, alpha=0.9)
    ax.axvline(0, linestyle="--", color="#333333", linewidth=1)
    ax.set_xlabel(r"Estimated origin shift  $\hat{\epsilon}_i$  (years)")
    ax.set_ylabel("Number of subjects")
    ax.set_title("Estimated origin shifts", fontweight="bold", fontsize=9)
    ax = axes[1]
    ax.scatter(s["A_obs"], s["S_hat"], s=11, alpha=0.35, color=UIUC_BLUE,
               edgecolors="none")
    lim = [0, float(np.nanmax([s["A_obs"].max(), s["S_hat"].max()])) * 1.02]
    ax.plot(lim, lim, "--", color="#999999", linewidth=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_aspect("equal")
    ax.set_xlabel("Recorded landmark age (years)")
    ax.set_ylabel("Calibrated landmark age (years)")
    ax.set_title("Recorded vs calibrated age", fontweight="bold", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig4a_shift_distribution.png"))
    fig.savefig(os.path.join(FIG, "fig4a_shift_distribution.pdf"))
    # standalone calibrated-age panel (backward-compatible filename)
    fig2, ax2 = plt.subplots(figsize=(4.6, 4.4))
    ax2.scatter(s["A_obs"], s["S_hat"], s=11, alpha=0.35, color=UIUC_BLUE,
                edgecolors="none")
    ax2.plot(lim, lim, "--", color="#999999", linewidth=1)
    ax2.set_xlim(lim); ax2.set_ylim(lim); ax2.set_aspect("equal")
    ax2.set_xlabel("Recorded landmark age (years)")
    ax2.set_ylabel("SILK-calibrated landmark age (years)")
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, "fig4b_calibrated_landmark.png"))
    fig2.savefig(os.path.join(FIG, "fig4b_calibrated_landmark.pdf"))
    plt.close(fig); plt.close(fig2)
    print("wrote fig4a / fig4b")


def fig5_km():
    """KM curves by SILK (Gaussian) 5-year risk tertile at the 3-year landmark."""
    pred = pd.read_csv(os.path.join(RES, "macs_predictions.csv"))
    z = pred[(pred["method"] == "SILK") & (pred["horizon"] == 5) &
             (pred["landmark_set"] == "fixed_3y")].copy()
    if len(z) < 30:
        z = pred[(pred["method"] == "SILK") & (pred["horizon"] == 5)].copy()
    z = z.drop_duplicates("subject_id")
    U = z["residual_time"].values.astype(float)
    D = z["event_status"].values.astype(int)
    p = z["risk_pred"].values.astype(float)
    ter = np.quantile(p, [1 / 3, 2 / 3])
    grp = np.where(p <= ter[0], 0, np.where(p <= ter[1], 1, 2))
    names = ["Low risk (bottom third)", "Medium risk", "High risk (top third)"]
    cols = [UIUC_BLUE, "#7F7F7F", UIUC_ORANGE]

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
    for g, nm, c in zip([0, 1, 2], names, cols):
        mask = grp == g
        t, s = km(U[mask], D[mask])
        ax.step(t, s, where="post", color=c, linewidth=2.0, label=nm)
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time from landmark (years)")
    ax.set_ylabel("AIDS-free survival probability")
    ax.legend(loc="lower left", frameon=False, fontsize=8.5)
    fig.tight_layout()
    fig.savefig(os.path.join(FIG, "fig5_km_risk_strata.png"))
    fig.savefig(os.path.join(FIG, "fig5_km_risk_strata.pdf"))
    plt.close(fig)
    print("wrote fig5_km_risk_strata")


if __name__ == "__main__":
    fig1_auc()
    fig2_brier()
    fig3_calibration()
    fig4_shifts()
    fig5_km()
    print("All MACS figures regenerated at %d dpi -> %s" % (DPI, FIG))
