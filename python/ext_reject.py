"""Why draws are rejected at the 2018 step (thesis Chapter 7): the 2018 level is either below the floor set by
non-social acquisition (prevalence exceeds 7% even with no social exposure) or above what any exposure rate up to
the solver bound can reach. Writes outputs/results_rejx.json (keys rejx_...)."""
from __future__ import annotations

import json, pathlib
import numpy as np

from dpm import calibrate as cal
from dpm.model import Params, equilibrium

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "outputs"


def main():
    rows = json.load(open(OUT / "calibration_draws.json"))
    low, high, g_low, g_high = 0, 0, [], []
    for r in rows:
        if r["ok"] or r["reason"] != "2018 level not reachable":
            continue
        base = Params(b=0.0, **{k: r[k] for k in cal.PARAM_ORDER})
        if equilibrium(base)[2] > r["t18"]:
            low += 1; g_low.append(r["gamma2"])
        elif equilibrium(base.with_(b=200.0))[2] < r["t18"]:
            high += 1; g_high.append(r["gamma2"])
    R = {"rejx_2018_floor_above": {"value": low, "note": "draws whose floor (b = 0) exceeds the 2018 target"},
         "rejx_2018_unreachable": {"value": high, "note": "draws whose equilibrium at b = 200 stays below the 2018 target"},
         "rejx_2018_unreachable_median_gamma2": {"value": float(np.median(g_high)) if g_high else float("nan"), "note": ""},
         "rejx_2018_unreachable_median_sigma": {"value": float(np.median([r["sigma"] for r in rows if not r["ok"] and r["reason"] == "2018 level not reachable"
                                                                          and equilibrium(Params(b=200.0, **{k: r[k] for k in cal.PARAM_ORDER}))[2] < r["t18"]])), "note": ""}}
    json.dump(R, open(OUT / "results_rejx.json", "w"), indent=1)
    print({k: v["value"] for k, v in R.items()})


if __name__ == "__main__":
    main()
