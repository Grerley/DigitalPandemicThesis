// agesex.js: Screen 8, the age- and sex-structured model (thesis Chapters 4, 5, 7 and 9).

import { h } from "../dom.js";
import { panel, metricGrid, chartCard, legend, condNote } from "../ui.js";
import { timeSeries } from "../charts.js";
import { viewHead, table, source } from "../common.js";
import { LAB, R } from "../../model/params.js";

export function AgeSexView() {
  const el = h("div", { class: "view" });
  el.appendChild(viewHead("Screen 8 · Age and sex", "Girls, boys and age bands",
    "A structured version splits young people into eight groups by sex and four age bands (10 to 12, 13 to 15, 16 to 19 and 20 to 24), with ageing through the bands, mixing mostly within age band, " +
    "and a vulnerability multiplier for girls and for the youngest band. It is fitted exactly to the 2022 HBSC figures for all young people, girls (13%), boys (9%) and 11-year-olds (9%), and then compared with what is known about 2018, which it was not fitted to."));

  const A = LAB.agesex, d18g = R("asxs_derived18_girls") / 100, d18b = R("asxs_derived18_boys") / 100;
  const mk = (which, title) => chartCard({ title, sub: "HBSC ages (11, 13 and 15), central fit (thesis results)", filename: `agesex-${which}`,
    render: () => timeSeries({ w: 560, h: 420, xs: A.years, yTitle: "Problematic use", yMax: 0.18,
      series: [{ key: "g", name: "Girls", color: "var(--c-i)", ys: A[which].girls.map((v) => v / 100) },
        { key: "b", name: "Boys", color: "var(--c-s)", dash: "6 4", ys: A[which].boys.map((v) => v / 100) }],
      points: [{ x: 2022, y: 0.13, label: "13%", color: "var(--c-i)" }, { x: 2022, y: 0.09, label: "9%", color: "var(--c-s)" },
        { x: 2018, y: d18g, label: "girls 2018*", color: "var(--c-i)" }, { x: 2018, y: d18b, label: "boys 2018*", color: "var(--c-s)" }] }) });

  el.appendChild(panel({ title: "Two accounts of the sex gap", children: [
    legend([{ name: "Girls", color: "var(--c-i)" }, { name: "Boys", color: "var(--c-s)", dash: "6 4" }]),
    h("div", { class: "grid cols-2" }, [mk("const", "Constant vulnerability of girls"), mk("shift", "Girls' vulnerability changes in 2018")]),
    h("p", { class: "hint", text: "* 2018 values derived from the pooled odds ratios of Boer et al. (2022), not published averages." }),
    metricGrid([
      { label: "Vulnerability of girls relative to boys", value: R("asx_d_P_phi_median").toFixed(2), foot: `median across draws (${R("asx_d_P_phi_p05").toFixed(2)} to ${R("asx_d_P_phi_p95").toFixed(2)})`, accent: true },
      { label: "2018 girls-to-boys ratio implied", value: R("p5_P_ratio18_sex_median").toFixed(2), foot: `constant vulnerability (${R("p5_P_ratio18_sex_p05").toFixed(2)} to ${R("p5_P_ratio18_sex_p95").toFixed(2)})` },
      { label: "Ratio from pooled odds ratios", value: R("asxs_derived18_ratio_sex").toFixed(2), foot: "2018 benchmark" },
      { label: "Italian national HBSC report", value: "about 1.7 to 1.8", foot: "2018, ages 13 and 15" },
    ]),
    h("p", { class: "prose", html: "Reproducing the 2022 breakdowns requires girls to be more vulnerable to acquisition than boys, almost whatever the assumed mixing. Tested against 2018, the model reproduces the age gradient. " +
      "On sex, a constant difference in vulnerability implies a 2018 ratio well above the one derived from pooled odds ratios, while its range overlaps the ratios in the Italian national report. " +
      "The verdict therefore depends on a benchmark that tabulations from the HBSC 2017/2018 open-access data could settle. If the pooled benchmark holds, girls' vulnerability rose after 2018, and girls account for about two thirds of the absolute rise." }),
    source("thesis Chapters 4, 5 and 7 and the age and sex paper (Mutibura and Lubuma, manuscript in preparation); Boniel-Nissim et al. (2024); Boer et al. (2022); Marino et al. (2020); python/ext_agesex.py, ext_agesex_shift.py and ext_paper5.py."),
  ] }));

  const bands = [["10 to 12", "b1012"], ["13 to 15", "b1315"], ["16 to 19", "b1619"], ["20 to 24", "b2024"], ["HBSC ages", "hbsc"], ["All 10 to 24", "whole"]];
  el.appendChild(panel({ title: "A minimum age of 16, represented by age band", sub: `relative reduction in problematic use, from ${R("asx_start_year")}`, children: [
    table(["Group", "Compliance 30%", "Compliance 50%", "Compliance 70%"],
      bands.map(([t, k]) => [t, ...[30, 50, 70].map((c) => `${R(`asx_ma_c${c}_${k}_red_pct`).toFixed(0)}%`)]), ["l", "r", "r", "r"]),
    h("p", { class: "prose", html: "Problematic use among the under-16s falls roughly in proportion to compliance, largely by construction. The effect fades among 16-to-24-year-olds. Under the model's assumption that rates after 16 do not depend on age at first access, delaying access delays problematic use more than it prevents it." }),
    condNote("Compliance levels are assumptions. South Africa has no enacted minimum age for social media; the scenario is hypothetical for South Africa."),
  ] }));
  return { el };
}
