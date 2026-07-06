// format.js — Number formatting for legible, tabular metrics.
// -----------------------------------------------------------------------------

export const pct = (x, dp = 1) => `${(x * 100).toFixed(dp)}%`;
export const pct0 = (x) => `${Math.round(x * 100)}%`;
export const fx = (x, dp = 2) => Number(x).toFixed(dp);
export const sig = (x, n = 3) => {
  if (x === 0) return "0";
  const d = Math.ceil(Math.log10(Math.abs(x)));
  const power = n - d;
  const mag = Math.pow(10, power);
  return String(Math.round(x * mag) / mag);
};
export const yr = (x, dp = 1) => `${x.toFixed(dp)} yr`;
export const signedPct = (x, dp = 1) => `${x >= 0 ? "−" : "+"}${Math.abs(x * 100).toFixed(dp)}%`;
export const zar = (x) => `R${x.toLocaleString("en-ZA", { maximumFractionDigits: 0 })}`;
