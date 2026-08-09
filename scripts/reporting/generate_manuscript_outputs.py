"""Generate manuscript tables and publication figures from existing result CSVs.

Reporting only: this script does not fit models or calculate new inferential
statistics.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


PROJECT_ROOT = Path(__file__).resolve().parents[2]
RESULT_ROOT = PROJECT_ROOT / "results"
TABLE_ROOT = PROJECT_ROOT / "results" / "tables"
FIGURE_ROOT = PROJECT_ROOT / "results" / "figures"
TABLE_ROOT.mkdir(parents=True, exist_ok=True)
FIGURE_ROOT.mkdir(parents=True, exist_ok=True)

f1 = pd.read_csv(RESULT_ROOT / "F1_weighted_descriptive.csv")
f2 = pd.read_csv(RESULT_ROOT / "F2_svy_logistic_frailty.csv")
f3 = pd.read_csv(RESULT_ROOT / "F3_reporting_scales.csv")
f4_logistic = pd.read_csv(RESULT_ROOT / "F4_sensitivity_logistic.csv")
f4_cox = pd.read_csv(RESULT_ROOT / "F4_sensitivity_cox.csv")
f4_no_fi = pd.read_csv(RESULT_ROOT / "F4_mortality_no_FI.csv")
f4_quartile_frailty = pd.read_csv(RESULT_ROOT / "F4_exposure_quartile_frailty.csv")
f4_quartile_mortality = pd.read_csv(RESULT_ROOT / "F4_exposure_quartile_mortality_no_FI.csv")
selection_comparison = pd.read_csv(RESULT_ROOT / "selection_comparison.csv")
visualization_quartiles = pd.read_csv(RESULT_ROOT / "result_visualization_quartiles.csv")
coverage = pd.read_csv(PROJECT_ROOT / "logs" / "frailty_component_coverage_audit.csv")
missingness = pd.read_csv(RESULT_ROOT / "analysis_missingness.csv")
contract = pd.read_csv(RESULT_ROOT / "analysis_contract.csv")
contract_values = dict(zip(contract.metric.astype(str), contract.value.astype(str)))


def contract_int(metric):
    return int(float(contract_values[metric]))


def format_ci(row, estimate="estimate", low="ci_low", high="ci_high"):
    return f"{row[estimate]:.3f} ({row[low]:.3f}-{row[high]:.3f})"


table1 = f1.copy()
table1["estimate_ci"] = table1.apply(
    lambda row: f"{row.estimate:.3f} ({row.ci_low:.3f}-{row.ci_high:.3f})", axis=1
)
table1.to_csv(TABLE_ROOT / "table1_weighted_descriptive.csv", index=False)

frailty_main = f2[f2.term == "diet_variety_5"].iloc[0]
mortality_main = f3[f3.term == "diet_variety_5"].iloc[0]
fi_mortality = f3[f3.term == "FI_primary"].iloc[0]

table2 = pd.DataFrame([
    {
        "analysis": "Frailty, primary survey-weighted logistic",
        "exposure": "DR1TNUMF per 5 items",
        "n": contract_int("primary_model_n"),
        "events": contract_int("primary_model_frail_events"),
        "effect": "OR",
        "estimate": frailty_main.odds_ratio,
        "ci_low": frailty_main.odds_ratio_ci_low,
        "ci_high": frailty_main.odds_ratio_ci_high,
        "p_value": frailty_main.p_value,
    },
    {
        "analysis": "Mortality, primary survey-weighted Cox",
        "exposure": "DR1TNUMF per 5 items",
        "n": contract_int("mortality_model_n"),
        "events": contract_int("mortality_model_events"),
        "effect": "HR",
        "estimate": mortality_main.hazard_ratio_reporting,
        "ci_low": mortality_main.hazard_ratio_ci_low_reporting,
        "ci_high": mortality_main.hazard_ratio_ci_high_reporting,
        "p_value": mortality_main.p_value,
    },
    {
        "analysis": "Mortality, primary survey-weighted Cox",
        "exposure": "FI per 0.1 unit",
        "n": contract_int("mortality_model_n"),
        "events": contract_int("mortality_model_events"),
        "effect": "HR",
        "estimate": fi_mortality.hazard_ratio_reporting,
        "ci_low": fi_mortality.hazard_ratio_ci_low_reporting,
        "ci_high": fi_mortality.hazard_ratio_ci_high_reporting,
        "p_value": fi_mortality.p_value,
    },
])
table2["effect_ci"] = table2.apply(format_ci, axis=1)
table2.to_csv(TABLE_ROOT / "table2_main_results.csv", index=False)

logistic = f4_logistic.copy()
logistic["analysis"] = logistic.model.map({
    "F4_same0_WTDRD1": "HUQ020 same=0",
    "F4_same1_WTDRD1": "HUQ020 same=1",
    "F4_primary_WTMEC2YR": "WTMEC2YR sensitivity",
    "F4_primary_energy_WTDRD1": "Total-energy adjustment",
})
logistic = logistic.drop(columns=["estimate"])
logistic = logistic.rename(columns={
    "odds_ratio": "estimate",
    "odds_ratio_ci_low": "ci_low",
    "odds_ratio_ci_high": "ci_high",
})
logistic = logistic.loc[:, ~logistic.columns.duplicated()]
logistic["effect"] = "OR"
logistic["effect_ci"] = logistic.apply(format_ci, axis=1)

cox = f4_cox[f4_cox.term == "diet_variety_5"].copy()
cox["analysis"] = cox.model.map({
    "F4_primary_PERMTH_INT": "Interview follow-up",
    "F4_primary_energy_WTDRD1": "Total-energy adjustment",
    "F4_mortality_parsimonious_PERMTH_EXM": "Parsimonious Cox",
})
cox = cox.drop(columns=["estimate"])
cox = cox.rename(columns={
    "hazard_ratio": "estimate",
    "hazard_ratio_ci_low": "ci_low",
    "hazard_ratio_ci_high": "ci_high",
})
cox = cox.loc[:, ~cox.columns.duplicated()]
cox["effect"] = "HR"
cox["effect_ci"] = cox.apply(format_ci, axis=1)

no_fi = f4_no_fi.copy()
no_fi["analysis"] = "No-FI Cox"
no_fi["effect"] = "HR"
no_fi = no_fi.drop(columns=["estimate"])
no_fi = no_fi.rename(columns={
    "hazard_ratio": "estimate",
    "hazard_ratio_ci_low": "ci_low",
    "hazard_ratio_ci_high": "ci_high",
})
no_fi["effect_ci"] = no_fi.apply(format_ci, axis=1)

quartile_frailty = f4_quartile_frailty.copy()
quartile_frailty["analysis"] = quartile_frailty.term.str.replace("exposure_q4", "Frailty ", regex=False) + " vs Q1"
quartile_frailty["effect"] = "OR"
quartile_frailty = quartile_frailty.rename(columns={
    "odds_ratio": "estimate",
    "odds_ratio_ci_low": "ci_low",
    "odds_ratio_ci_high": "ci_high",
})
quartile_frailty["effect_ci"] = quartile_frailty.apply(format_ci, axis=1)

quartile_mortality = f4_quartile_mortality.copy()
quartile_mortality["analysis"] = quartile_mortality.term.str.replace("exposure_q4", "Mortality ", regex=False) + " vs Q1"
quartile_mortality["effect"] = "HR"
quartile_mortality = quartile_mortality.rename(columns={
    "hazard_ratio": "estimate",
    "hazard_ratio_ci_low": "ci_low",
    "hazard_ratio_ci_high": "ci_high",
})
quartile_mortality["effect_ci"] = quartile_mortality.apply(format_ci, axis=1)

table3 = pd.concat([
    logistic[["analysis", "effect", "estimate", "ci_low", "ci_high", "p_value", "n", "events", "effect_ci"]],
    cox[["analysis", "effect", "estimate", "ci_low", "ci_high", "p_value", "n", "events", "effect_ci"]],
    no_fi[["analysis", "effect", "estimate", "ci_low", "ci_high", "p_value", "n", "events", "effect_ci"]],
    quartile_frailty[["analysis", "effect", "estimate", "ci_low", "ci_high", "p_value", "n", "events", "effect_ci"]],
    quartile_mortality[["analysis", "effect", "estimate", "ci_low", "ci_high", "p_value", "n", "events", "effect_ci"]],
], ignore_index=True)
table3.to_csv(TABLE_ROOT / "table3_sensitivity_results.csv", index=False)
coverage.to_csv(TABLE_ROOT / "tableS1_fi_component_coverage.csv", index=False)
missingness.to_csv(TABLE_ROOT / "tableS2_analysis_missingness.csv", index=False)
selection_comparison.to_csv(TABLE_ROOT / "tableS3_selection_comparison.csv", index=False)
pd.concat([
    f4_quartile_frailty.assign(outcome_label="frailty"),
    f4_quartile_mortality.assign(outcome_label="mortality_no_FI"),
], ignore_index=True).to_csv(TABLE_ROOT / "tableS4_exposure_quartile_results.csv", index=False)


plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.8,
    "savefig.facecolor": "white",
})


def interval_plot(frame, colors, xlabel, output, xlim, style="forest"):
    frame = frame.reset_index(drop=True)
    y = np.arange(len(frame))[::-1]
    fig, ax = plt.subplots(figsize=(7.1, 3.7))
    for index, row in frame.iterrows():
        color = colors[index]
        ypos = y[index]
        ax.plot([row.ci_low, row.ci_high], [ypos, ypos], color=color, linewidth=1.7, solid_capstyle="round")
        ax.plot(row.estimate, ypos, "o", color=color, markersize=6.5, markeredgecolor="white", markeredgewidth=0.8)
    ax.axvline(1.0, color="#777777", linewidth=0.85, linestyle=(0, (3, 3)), zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(*xlim)
    ax.set_yticks(y)
    ax.set_yticklabels(frame.label)
    ax.set_xlabel(xlabel)
    ax.grid(axis="x", color="#d9d9d9", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#555555")
    ax.spines["bottom"].set_color("#555555")
    if style == "grouped":
        for ypos in y[1::2]:
            ax.axhspan(ypos - 0.42, ypos + 0.42, color="#f5f5f2", zorder=-1)
    fig.tight_layout(pad=0.8)
    fig.savefig(output, dpi=400, bbox_inches="tight")
    plt.close(fig)


# Figure-label helper is defined before the multi-panel result figures.
def add_panel_label(ax, label, y=0.98):
    ax.text(0.02, y, label, transform=ax.transAxes, fontsize=12,
            fontweight="bold", va="top", ha="left", color="#202020",
            bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5}, zorder=5)


# Figure 1: primary model results. This is the main inferential result figure,
# not a collection of sensitivity specifications.
primary_effects = pd.DataFrame([
    {
        "label": "Frailty per 5 items",
        "effect": "OR",
        "estimate": frailty_main.odds_ratio,
        "ci_low": frailty_main.odds_ratio_ci_low,
        "ci_high": frailty_main.odds_ratio_ci_high,
        "color": "#173f5f",
    },
    {
        "label": "Mortality per 5 items",
        "effect": "HR",
        "estimate": mortality_main.hazard_ratio_reporting,
        "ci_low": mortality_main.hazard_ratio_ci_low_reporting,
        "ci_high": mortality_main.hazard_ratio_ci_high_reporting,
        "color": "#b84a62",
    },
    {
        "label": "Mortality per 0.1 FI",
        "effect": "HR",
        "estimate": fi_mortality.hazard_ratio_reporting,
        "ci_low": fi_mortality.hazard_ratio_ci_low_reporting,
        "ci_high": fi_mortality.hazard_ratio_ci_high_reporting,
        "color": "#9a7b35",
    },
])
fig, ax = plt.subplots(figsize=(7.1, 3.15))
y = np.arange(len(primary_effects))[::-1]
for idx, row in primary_effects.reset_index(drop=True).iterrows():
    ypos = y[idx]
    ax.plot([row.ci_low, row.ci_high], [ypos, ypos], color=row.color,
            linewidth=2.2, solid_capstyle="round")
    ax.plot(row.estimate, ypos, "o", color=row.color, markersize=8,
            markeredgecolor="white", markeredgewidth=0.9)
    ax.text(2.35, ypos, f"{row.effect} {row.estimate:.3f} ({row.ci_low:.3f}-{row.ci_high:.3f})",
            va="center", ha="left", fontsize=8.7, color="#38454d")
ax.axvline(1.0, color="#777777", linewidth=0.9, linestyle=(0, (3, 3)))
ax.set_xscale("log")
ax.set_xlim(0.65, 3.2)
ax.set_yticks(y)
ax.set_yticklabels(primary_effects.label)
ax.set_xlabel("Effect estimate (log scale; 95% CI)")
ax.set_title("Primary model results", fontsize=12, fontweight="bold", color="#173f5f", pad=8)
ax.grid(axis="x", color="#d9d9d9", linewidth=0.55)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#555555")
ax.spines["bottom"].set_color("#555555")
fig.tight_layout(pad=0.8)
fig.savefig(FIGURE_ROOT / "figure1_primary_results.png", dpi=400, bbox_inches="tight")
plt.close(fig)


# Figure 2: weighted outcome proportions across exposure quartiles. Both
# panels are descriptive outcome results and retain the survey-weighted CIs.
quartile_plot = visualization_quartiles.copy()
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.55), sharey=False)
panel_specs = [
    ("frailty_prevalence", "Frailty prevalence", "Weighted proportion", "#173f5f", "A"),
    ("linked_mortality_proportion", "Linked mortality", "Weighted proportion", "#b84a62", "B"),
]
for ax, (outcome, title, ylabel, color, panel) in zip(axes, panel_specs):
    frame = quartile_plot[quartile_plot.outcome == outcome].copy()
    frame["exposure_q4"] = pd.Categorical(frame.exposure_q4, categories=["Q1", "Q2", "Q3", "Q4"], ordered=True)
    frame = frame.sort_values("exposure_q4")
    x = np.arange(len(frame))
    ax.bar(x, frame.weighted_estimate, color=color, alpha=0.82, width=0.66,
           edgecolor="white", linewidth=0.7)
    ax.errorbar(x, frame.weighted_estimate,
                yerr=[frame.weighted_estimate - frame.ci_low, frame.ci_high - frame.weighted_estimate],
                fmt="none", ecolor="#202020", elinewidth=1.0, capsize=3)
    for xpos, (_, row) in zip(x, frame.iterrows()):
        ax.text(xpos, min(0.99, row.ci_high + 0.035), f"n={int(row.n)}\ne={int(row.events)}",
                ha="center", va="bottom", fontsize=7.8, color="#38454d")
    ax.set_xticks(x, ["Q1", "Q2", "Q3", "Q4"])
    ax.set_ylim(0, min(1.0, max(0.5, float(frame.ci_high.max()) + 0.20)))
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10.2, pad=7)
    ax.grid(axis="y", color="#dedede", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    add_panel_label(ax, panel)
fig.suptitle("Weighted outcome proportions by exposure quartile", fontsize=11.5,
             fontweight="bold", y=1.02, color="#173f5f")
fig.tight_layout(pad=0.75)
fig.savefig(FIGURE_ROOT / "figure2_quartile_outcomes.png", dpi=400, bbox_inches="tight")
plt.close(fig)


frailty_plot = pd.DataFrame([
    {"label": "Primary WTDRD1", "estimate": frailty_main.odds_ratio, "ci_low": frailty_main.odds_ratio_ci_low, "ci_high": frailty_main.odds_ratio_ci_high},
    {"label": "HUQ020 same=0", "estimate": logistic.loc[logistic.analysis == "HUQ020 same=0", "estimate"].iloc[0], "ci_low": logistic.loc[logistic.analysis == "HUQ020 same=0", "ci_low"].iloc[0], "ci_high": logistic.loc[logistic.analysis == "HUQ020 same=0", "ci_high"].iloc[0]},
    {"label": "HUQ020 same=1", "estimate": logistic.loc[logistic.analysis == "HUQ020 same=1", "estimate"].iloc[0], "ci_low": logistic.loc[logistic.analysis == "HUQ020 same=1", "ci_low"].iloc[0], "ci_high": logistic.loc[logistic.analysis == "HUQ020 same=1", "ci_high"].iloc[0]},
    {"label": "Energy-adjusted", "estimate": logistic.loc[logistic.analysis == "Total-energy adjustment", "estimate"].iloc[0], "ci_low": logistic.loc[logistic.analysis == "Total-energy adjustment", "ci_low"].iloc[0], "ci_high": logistic.loc[logistic.analysis == "Total-energy adjustment", "ci_high"].iloc[0]},
    {"label": "WTMEC2YR", "estimate": logistic.loc[logistic.analysis == "WTMEC2YR sensitivity", "estimate"].iloc[0], "ci_low": logistic.loc[logistic.analysis == "WTMEC2YR sensitivity", "ci_low"].iloc[0], "ci_high": logistic.loc[logistic.analysis == "WTMEC2YR sensitivity", "ci_high"].iloc[0]},
])
interval_plot(
frailty_plot,
    ["#173f5f", "#4f6d7a", "#4f6d7a", "#b84a62", "#9a7b35"],
    "Odds ratio (log scale)",
    FIGURE_ROOT / "figure6_frailty_sensitivity.png",
    (0.65, 1.25),
    style="forest",
)

mortality_plot = cox[["analysis", "estimate", "ci_low", "ci_high"]].rename(columns={"analysis": "label"})
primary_mortality = pd.DataFrame([{
    "label": "Primary PERMTH_EXM",
    "estimate": mortality_main.hazard_ratio_reporting,
    "ci_low": mortality_main.hazard_ratio_ci_low_reporting,
    "ci_high": mortality_main.hazard_ratio_ci_high_reporting,
}])
mortality_plot = pd.concat([primary_mortality, mortality_plot, no_fi[["analysis", "estimate", "ci_low", "ci_high"]].rename(columns={"analysis": "label"})], ignore_index=True)
interval_plot(
    mortality_plot,
    ["#173f5f", "#7a5195", "#b84a62", "#3f7d5a", "#5b6770"],
    "Hazard ratio (log scale)",
    FIGURE_ROOT / "figure7_mortality_sensitivity.png",
    (0.6, 1.8),
    style="grouped",
)


# Figure 3: weighted profiles across exposure quartiles. The dot-and-whisker
# layout is intentionally distinct from the connected trajectory in Figure 5.
profile_specs = [
    ("fi_mean", "Frailty Index", "Weighted mean FI", "#b84a62", "A"),
    ("exposure_mean", "Reported food/beverage count", "Weighted mean count", "#4f6d7a", "B"),
]
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.55), sharey=False)
for ax, (outcome, title, ylabel, color, panel) in zip(axes, profile_specs):
    frame = quartile_plot[quartile_plot.outcome == outcome].copy()
    frame["exposure_q4"] = pd.Categorical(frame.exposure_q4, categories=["Q1", "Q2", "Q3", "Q4"], ordered=True)
    frame = frame.sort_values("exposure_q4")
    y = np.arange(len(frame))[::-1]
    ax.errorbar(frame.weighted_estimate, y,
                xerr=[frame.weighted_estimate - frame.ci_low, frame.ci_high - frame.weighted_estimate],
                fmt="o", color=color, ecolor="#202020", markersize=6.5,
                markeredgecolor="white", markeredgewidth=0.8, capsize=3, linewidth=1.1)
    ax.set_yticks(y, ["Q1", "Q2", "Q3", "Q4"])
    ax.set_xlabel(ylabel)
    ax.set_title(title, fontsize=10.0, pad=8)
    ax.grid(axis="x", color="#dedede", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    add_panel_label(ax, panel)
fig.suptitle("Weighted exposure and FI profiles by quartile", fontsize=11.5,
             fontweight="bold", y=1.02, color="#173f5f")
fig.tight_layout(pad=0.75)
fig.savefig(FIGURE_ROOT / "figure3_weighted_profiles.png", dpi=400, bbox_inches="tight")
plt.close(fig)


# Figure 4: selection audit. Raw scales are kept in separate panels because
# the variables have different units; missingness is shown for the excluded group.
selection_numeric = selection_comparison[
    selection_comparison.variable.isin(["RIDAGEYR", "DR1TNUMF", "FI_primary", "INDFMPIR"])
].copy()
metric_specs = [
    ("RIDAGEYR", "Age (years)", "#173f5f"),
    ("DR1TNUMF", "Reported food/beverage count", "#4f6d7a"),
    ("FI_primary", "Frailty Index", "#b84a62"),
    ("INDFMPIR", "Poverty-income ratio", "#9a7b35"),
]
fig, axes = plt.subplots(2, 3, figsize=(7.1, 5.1))
axes = axes.ravel()
for idx, (variable, title, color) in enumerate(metric_specs):
    ax = axes[idx]
    rows = selection_numeric[selection_numeric.variable == variable].set_index("group")
    complete = rows.loc["primary_complete_case", "mean_or_proportion"]
    excluded = rows.loc["excluded_from_primary", "mean_or_proportion"]
    complete_sd = rows.loc["primary_complete_case", "standard_deviation"]
    excluded_sd = rows.loc["excluded_from_primary", "standard_deviation"]
    ax.plot([0, 1], [complete, excluded], color="#a8afb3", linewidth=1.4, zorder=1)
    ax.errorbar([0], [complete], yerr=[complete_sd], fmt="o", color=color,
                markersize=6.5, capsize=3, linewidth=1.1, zorder=2, label="Complete")
    ax.errorbar([1], [excluded], yerr=[excluded_sd], fmt="o", color="#202020",
                markersize=6.5, capsize=3, linewidth=1.1, zorder=2, label="Excluded")
    ax.set_title(title, fontsize=9.6, pad=7)
    ax.set_xticks([0, 1], ["Complete", "Excluded"])
    ax.grid(axis="y", color="#dedede", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    add_panel_label(ax, chr(ord("A") + idx))

ax = axes[4]
missing_rows = selection_comparison[
    (selection_comparison.group == "excluded_from_primary") &
    selection_comparison.variable.isin(["DR1TNUMF", "INDFMPIR", "FI_primary", "WTDRD1"])
].copy()
missing_rows["label"] = missing_rows.variable.map({
    "DR1TNUMF": "Food/beverage count",
    "INDFMPIR": "Poverty-income ratio",
    "FI_primary": "Frailty Index",
    "WTDRD1": "Dietary-day weight",
})
missing_rows = missing_rows.set_index("label").reindex([
    "Food/beverage count", "Poverty-income ratio", "Frailty Index", "Dietary-day weight"
]).reset_index()
ax.barh(missing_rows.label, missing_rows.missing_percent, color="#b84a62", alpha=0.82)
ax.set_xlabel("Missing among excluded (%)")
ax.invert_yaxis()
ax.grid(axis="x", color="#dedede", linewidth=0.55)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#777777")
ax.spines["bottom"].set_color("#777777")
add_panel_label(ax, "E")

ax = axes[5]
sex_rows = selection_comparison[
    (selection_comparison.variable == "RIAGENDR") &
    selection_comparison.level.isin(["1", "2"])
].copy()
sex_pivot = sex_rows.pivot(index="group", columns="level", values="mean_or_proportion").reindex([
    "primary_complete_case", "excluded_from_primary"
])
bottom = np.zeros(2)
sex_colors = {"1": "#173f5f", "2": "#d7a33d"}
sex_labels = {"1": "Male", "2": "Female"}
for level in ["1", "2"]:
    values = sex_pivot[level].to_numpy(dtype=float)
    ax.bar([0, 1], values, bottom=bottom, color=sex_colors[level], width=0.62,
           edgecolor="white", linewidth=0.7, label=sex_labels[level])
    for xpos, value, base in zip([0, 1], values, bottom):
        ax.text(xpos, base + value / 2, f"{value * 100:.0f}%", ha="center", va="center",
                fontsize=8.2, color="white" if level == "1" else "#202020")
    bottom += values
ax.set_xticks([0, 1], ["Keep", "Drop"])
ax.set_ylim(0, 1)
ax.set_ylabel("Sex composition")
ax.yaxis.set_major_formatter(lambda value, position: f"{value * 100:.0f}%")
ax.set_title("Sex composition", fontsize=9.6, pad=7)
ax.legend(frameon=False, fontsize=7.2, loc="upper center", bbox_to_anchor=(0.5, -0.18), ncol=2)
ax.grid(axis="y", color="#dedede", linewidth=0.55)
ax.set_axisbelow(True)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_color("#777777")
ax.spines["bottom"].set_color("#777777")
add_panel_label(ax, "F")
fig.suptitle("Complete-case selection and missingness audit", fontsize=12, fontweight="bold", y=1.03, color="#173f5f")
fig.tight_layout(pad=0.9, rect=[0, 0.08, 1, 0.96])
fig.savefig(FIGURE_ROOT / "figure4_selection_audit.png", dpi=400, bbox_inches="tight")
plt.close(fig)


# Figure 5: quartile trajectories. The connected estimates make the ordered
# frailty pattern and the non-monotonic mortality pattern directly comparable.
frailty_quartile = {"Q1": (1.0, 1.0, 1.0)}
for _, row in f4_quartile_frailty.iterrows():
    q = row.term[-2:]
    frailty_quartile[q] = (row.odds_ratio, row.odds_ratio_ci_low, row.odds_ratio_ci_high)
mortality_quartile = {"Q1": (1.0, 1.0, 1.0)}
for _, row in f4_quartile_mortality.iterrows():
    q = row.term[-2:]
    mortality_quartile[q] = (row.hazard_ratio, row.hazard_ratio_ci_low, row.hazard_ratio_ci_high)

fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.55), sharex=True)
for ax, values, ylabel, title, color, panel in [
    (axes[0], frailty_quartile, "Odds ratio (log scale)", "Frailty", "#173f5f", "A"),
    (axes[1], mortality_quartile, "Hazard ratio (log scale)", "Mortality without FI", "#b84a62", "B"),
]:
    x = np.arange(1, 5)
    est = np.array([values[f"Q{i}"][0] for i in range(1, 5)])
    low = np.array([values[f"Q{i}"][1] for i in range(1, 5)])
    high = np.array([values[f"Q{i}"][2] for i in range(1, 5)])
    # Q1 is the reference category and has no estimated confidence interval.
    ax.fill_between(x[1:], low[1:], high[1:], color=color, alpha=0.14, linewidth=0)
    ax.plot(x, est, color=color, linewidth=2.0, marker="o", markersize=5.6,
            markeredgecolor="white", markeredgewidth=0.8)
    ax.axhline(1.0, color="#777777", linewidth=0.8, linestyle=(0, (3, 3)))
    ax.set_yscale("log")
    ax.set_xticks(x, ["Q1", "Q2", "Q3", "Q4"])
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=10.2, pad=7)
    ax.grid(axis="y", color="#dedede", linewidth=0.55)
    ax.set_axisbelow(True)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color("#777777")
    ax.spines["bottom"].set_color("#777777")
    add_panel_label(ax, panel)
axes[0].text(1, 0.31, "Reference", fontsize=8.2, color="#5b6770", ha="center")
axes[1].text(1, 0.31, "Reference", fontsize=8.2, color="#5b6770", ha="center")
fig.suptitle("Exposure-quartile sensitivity: effect estimates versus Q1", fontsize=11.5,
             fontweight="bold", y=1.02, color="#173f5f")
fig.tight_layout(pad=0.7)
fig.savefig(FIGURE_ROOT / "figure5_quartile_trajectories.png", dpi=400, bbox_inches="tight")
plt.close(fig)

manifest = pd.DataFrame([
    {"artifact": "table1_weighted_descriptive.csv", "source": "F1_weighted_descriptive.csv"},
    {"artifact": "table2_main_results.csv", "source": "F2_svy_logistic_frailty.csv; F3_reporting_scales.csv"},
    {"artifact": "table3_sensitivity_results.csv", "source": "F4_sensitivity_logistic.csv; F4_sensitivity_cox.csv; F4_mortality_no_FI.csv"},
    {"artifact": "tableS1_fi_component_coverage.csv", "source": "frailty_component_coverage_audit.csv"},
    {"artifact": "tableS2_analysis_missingness.csv", "source": "analysis_missingness.csv"},
    {"artifact": "tableS3_selection_comparison.csv", "source": "selection_comparison.csv"},
    {"artifact": "tableS4_exposure_quartile_results.csv", "source": "F4_exposure_quartile_frailty.csv; F4_exposure_quartile_mortality_no_FI.csv"},
    {"artifact": "figure1_primary_results.png", "source": "F2_svy_logistic_frailty.csv; F3_reporting_scales.csv"},
    {"artifact": "figure2_quartile_outcomes.png", "source": "result_visualization_quartiles.csv"},
    {"artifact": "figure3_weighted_profiles.png", "source": "result_visualization_quartiles.csv"},
    {"artifact": "figure4_selection_audit.png", "source": "selection_comparison.csv"},
    {"artifact": "figure5_quartile_trajectories.png", "source": "F4_exposure_quartile_frailty.csv; F4_exposure_quartile_mortality_no_FI.csv"},
    {"artifact": "figure6_frailty_sensitivity.png", "source": "F2_svy_logistic_frailty.csv; F4_sensitivity_logistic.csv"},
    {"artifact": "figure7_mortality_sensitivity.png", "source": "F3_reporting_scales.csv; F4_sensitivity_cox.csv; F4_mortality_no_FI.csv"},
])
manifest.to_csv(TABLE_ROOT / "manuscript_output_manifest.csv", index=False)
print("MANUSCRIPT_OUTPUTS_PASS")
print(f"tables={len(list(TABLE_ROOT.glob('*.csv')))}")
print(f"figures={len(list(FIGURE_ROOT.glob('*.png')))}")


