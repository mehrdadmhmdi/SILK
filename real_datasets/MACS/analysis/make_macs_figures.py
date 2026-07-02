#!/usr/bin/env python3
"""
Publication-quality MACS figures (600 dpi), rebuilt from the result CSVs.
Replaces the cramped ggplot output: readable fonts, distinguishable series,
non-truncated legend. Matches the paper's matplotlib figure style.

Usage:
    python3 make_macs_figures.py
Reads:  results/macs_metrics_per_horizon.csv, macs_predictions.csv,
        macs_silk_shifts.csv
Writes: figures/fig1..fig5  (.png and .pdf, 600 dpi)
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
DPI = 600

# ---- house style ------------------------------------------------------------
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": DPI,
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 11,
    "axes.edgecolor": "#444444",
    "axes.linewidth": 0.8,
    "axes.grid": True,
    "grid.color": "#DDDDDD",
    "grid.linewidth": 0.6,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "savefig.bbox": "tight",
})

UIUC_ORANGE = "#FF5F05"
UIUC_BLUE = "#13294B"

# display label -> (color, linestyle, marker)
STYLE = {
    "SILK Gaussian Kernel + Cox": (UIUC_ORANGE, "-",  "o"),
    "SILK Linear Kernel + Cox":   ("#8C564B",   "--", "s"),
    "Mixed-Model Landmarking":    (UIUC_BLUE,   "-.", "^"),
    "Joint Model":                ("#C44E52",   ":",  "D"),
    "Landmarking":                ("#7F7F7F",   "-",  "x"),
}
ORDER = list(STYLE.keys())
LSMAP = {  # raw method id -> display label
    "SILK": "SILK Gaussian Kernel + Cox",
    "SILK-LinearMMD": "SILK Linear Kernel + Cox",
    "MMLM-Recorded": "Mixed-Model Landmarking",
    "JM-Recorded": "Joint Model",
    "Landmark-Recorded": "Landmarking",
}
FACET_TITLE = {"fixed_1y": "1-year landmark", "fixed_2y": "2-year landmark",
               "fixed_3y": "3-year landmark", "fixed_5y": "5-year landmark"}


def _load_metrics():
    m = pd.read_csv(os.path.join(RES, "macs_metrics_per_horizon.csv"))
    m["label"] = m["method"].map(LSMAP)
    return m


def _facet_grid(metric, se_col, ylab, title, outfile, ylim=None):
    m = _load_metrics()
    sets = ["fixed_1y", "fixed_2y", "fixed_3y", "fixed_5y"]
    horizons = sorted(m["horizon"].unique())
    fig, axes = plt.subplots(2, 2, figsize=(9.5, 7.2), sharex=True, sharey=True)
    axes = axes.ravel()
    for ax, ls in zip(axes, sets):
        sub = m[m["landmark_set"] == ls]
        for lab in ORDER:
            d = sub[sub["label"] == lab].sort_values("horizon")
            if d.empty:
                continue
            color, style, marker = STYLE[lab]
            ax.plot(d["horizon"], d[metric], style, color=color, marker=marker,
                    markersize=6, linewidth=1.8, markerfacecolor=color if marker != "x" else color,
                    markeredgecolor=color, label=lab, zorder=3, clip_on=False)
            if se_col and se_col in d:
                ax.errorbar(d["horizon"], d[metric.replace("pooled", "mean")],
                            yerr=1.96 * d[se_col], fmt="none", ecolor=color,
                            elinewidth=0.9, capsize=2.5, alpha=0.55, zorder=2)
        ax.set_title(FACET_TITLE[ls], fontweight="bold")
        ax.set_xticks(horizons)
        if ylim:
            ax.set_ylim(*ylim)
        ax.margins(x=0.06)
    for ax in axes[2:]:
        ax.set_xlabel("Prediction horizon (years)")
    for ax in (axes[0], axes[2]):
        ax.set_ylabel(ylab)
    handles = [Line2D([0], [0], color=STYLE[l][0], linestyle=STYLE[l][1],
                      marker=STYLE[l][2], markersize=6, linewidth=1.8) for l in ORDER]
    fig.legend(handles, ORDER, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.tight_layout(rect=[0, 0.07, 1, 0.96])
    fig.savefig(os.path.join(FIG, outfile + ".png"))
    fig.savefig(os.path.join(FIG, outfile + ".pdf"))
    plt.close(fig)
    print("wrote", outfile)


def fig1_auc():
    _facet_grid("auc_pooled", "auc_se", "IPCW AUC",
                "Discrimination: time-dependent AUC by prediction horizon",
                "fig1_auc_by_horizon", ylim=(0.40, 0.92))


def fig2_brier():
    _facet_grid("brier_pooled", "brier_se", "IPCW Brier score",
                "Accuracy: IPCW Brier score by prediction horizon",
                "fig2_brier_by_horizon")


def fig3_calibration():
    """Predicted vs observed risk (deciles), pooled, at horizons 2 and 5 yr."""
    pred = pd.read_csv(os.path.join(RES, "macs_predictions.csv"))
    pred["at_risk"] = pred["at_risk"].astype(str).str.lower().isin(["true", "1"])
    # IPCW using KM of censoring on residual time
    from math import inf
    def km_censor(times, u, d):
        # survfit of censoring: event = 1-delta
        order = np.argsort(u)
        u_s = u[order]; c = (1 - d)[order]
        n = len(u_s); surv = 1.0; S = {}
        uniq = np.unique(u_s)
        at_risk = n; km = 1.0
        surv_at = {}
        # standard KM
        for t in uniq:
            dt = np.sum((u_s == t) & (c == 1))
            nt = np.sum(u_s >= t)
            if nt > 0:
                km *= (1 - dt / nt)
            surv_at[t] = km
        keys = np.array(sorted(surv_at))
        vals = np.array([surv_at[k] for k in keys])
        def G(x):
            out = np.empty(len(x))
            for i, xi in enumerate(x):
                idx = keys[keys <= xi]
                out[i] = vals[np.searchsorted(keys, idx[-1])] if len(idx) else 1.0
            return np.clip(out, 1e-6, 1)
        return G
    horizons = [2, 5]
    fig, axes = plt.subplots(1, 2, figsize=(9.5, 5.0))
    outc = pred.drop_duplicates("subject_id")[["subject_id", "residual_time", "event_status"]]
    Gfun = km_censor(outc["residual_time"].values.astype(float),
                     outc["residual_time"].values.astype(float),
                     outc["event_status"].values.astype(float))
    for ax, h in zip(axes, horizons):
        for lab in ORDER:
            raw = [k for k, v in LSMAP.items() if v == lab][0]
            z = pred[(pred["method"] == raw) & (pred["horizon"] == h)].copy()
            if len(z) < 40:
                continue
            U = z["residual_time"].values.astype(float)
            D = z["event_status"].values.astype(int)
            y = np.where((D == 1) & (U <= h), 1.0, np.where(U > h, 0.0, np.nan))
            w = np.where((D == 1) & (U <= h), 1 / Gfun(U),
                         np.where(U > h, 1 / Gfun(np.full(len(U), h)), 0.0))
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
            ax.plot(xs, ys, style, color=color, marker=marker, markersize=6,
                    linewidth=1.6, label=lab, zorder=3)
        amax = 0.6 if h == 5 else 0.35
        ax.plot([0, amax], [0, amax], "--", color="#888888", linewidth=1, zorder=1)
        ax.set_xlim(0, amax); ax.set_ylim(0, amax)
        ax.set_aspect("equal")
        ax.set_title(f"{h}-year horizon", fontweight="bold")
        ax.set_xlabel("Predicted risk")
    axes[0].set_ylabel("Observed proportion (IPCW)")
    handles = [Line2D([0], [0], color=STYLE[l][0], linestyle=STYLE[l][1],
                      marker=STYLE[l][2], markersize=6, linewidth=1.6) for l in ORDER]
    fig.legend(handles, ORDER, loc="lower center", ncol=3, frameon=False,
               bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("Calibration: predicted vs observed AIDS risk",
                 fontsize=13, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0.10, 1, 0.95])
    fig.savefig(os.path.join(FIG, "fig3_calibration.png"))
    fig.savefig(os.path.join(FIG, "fig3_calibration.pdf"))
    plt.close(fig)
    print("wrote fig3_calibration")


def fig4_shifts():
    s = pd.read_csv(os.path.join(RES, "macs_silk_shifts.csv"))
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.6))
    ax = axes[0]
    ax.hist(s["e_hat"], bins=40, color=UIUC_ORANGE, edgecolor="white", alpha=0.9)
    ax.axvline(0, linestyle="--", color="#333333", linewidth=1)
    ax.set_xlabel(r"Estimated origin shift  $\hat{\epsilon}_i$  (years)")
    ax.set_ylabel("Number of subjects")
    ax.set_title("SILK-estimated origin shifts", fontweight="bold")
    ax = axes[1]
    ax.scatter(s["A_obs"], s["S_hat"], s=14, alpha=0.35, color=UIUC_BLUE,
               edgecolors="none")
    lim = [0, np.nanmax([s["A_obs"].max(), s["S_hat"].max()]) * 1.02]
    ax.plot(lim, lim, "--", color="#888888", linewidth=1)
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Recorded landmark age (years since SC midpoint)")
    ax.set_ylabel("SILK-calibrated landmark age (years)")
    ax.set_title("Recorded vs calibrated landmark age", fontweight="bold")
    fig.suptitle("SILK registration on the MACS cohort (Gaussian kernel, full data)",
                 fontsize=13, fontweight="bold", y=1.0)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(os.path.join(FIG, "fig4a_shift_distribution.png"))
    fig.savefig(os.path.join(FIG, "fig4a_shift_distribution.pdf"))
    # also a standalone calibrated-age panel for backward-compatible filename
    fig2, ax2 = plt.subplots(figsize=(5.6, 5.2))
    ax2.scatter(s["A_obs"], s["S_hat"], s=14, alpha=0.35, color=UIUC_BLUE, edgecolors="none")
    ax2.plot(lim, lim, "--", color="#888888", linewidth=1)
    ax2.set_xlim(lim); ax2.set_ylim(lim)
    ax2.set_xlabel("Recorded landmark age (years)")
    ax2.set_ylabel("SILK-calibrated landmark age (years)")
    ax2.set_title("Recorded vs SILK-calibrated landmark age", fontweight="bold")
    fig2.tight_layout()
    fig2.savefig(os.path.join(FIG, "fig4b_calibrated_landmark.png"))
    fig2.savefig(os.path.join(FIG, "fig4b_calibrated_landmark.pdf"))
    plt.close(fig); plt.close(fig2)
    print("wrote fig4a / fig4b")


def fig5_km():
    """KM curves by SILK (Gaussian) risk tertile at the 5-year horizon,
    3-year landmark (largest well-populated risk set)."""
    pred = pd.read_csv(os.path.join(RES, "macs_predictions.csv"))
    z = pred[(pred["method"] == "SILK") & (pred["horizon"] == 5) &
             (pred["landmark_set"] == "fixed_3y")].copy()
    if len(z) < 30:
        z = pred[(pred["method"] == "SILK") & (pred["horizon"] == 5)].copy()
    z = z.drop_duplicates("subject_id")
    U = z["residual_time"].values.astype(float)
    D = z["event_status"].values.astype(int)
    p = z["risk_pred"].values.astype(float)
    ter = np.quantile(p, [1/3, 2/3])
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

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    xmax = min(10, np.nanmax(U))
    for g, nm, c in zip([0, 1, 2], names, cols):
        mask = grp == g
        t, s = km(U[mask], D[mask])
        ax.step(t, s, where="post", color=c, linewidth=2.2, label=nm)
    ax.set_xlim(0, xmax); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Time from landmark (years)")
    ax.set_ylabel("AIDS-free survival probability")
    ax.set_title("Kaplan–Meier survival by SILK risk tertile\n(3-year landmark)",
                 fontweight="bold")
    ax.legend(loc="lower left", frameon=False)
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
