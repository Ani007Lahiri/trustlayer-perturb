import anndata as ad, numpy as np
# Reads interval predictions + ground-truth solution, computes empirical coverage
# of the nominal-0.90 interval and the undercoverage gap. Metric is model-agnostic:
# it scores the UNCERTAINTY quality of any method that emits intervals.
par = {"input_prediction": None, "input_solution": None, "output": None}  # filled by viash
pred = ad.read_h5ad(par["input_prediction"])
sol  = ad.read_h5ad(par["input_solution"])
lo = pred.layers["interval_lower"]; hi = pred.layers["interval_upper"]
y  = sol.X
inside = (y >= lo) & (y <= hi)
coverage = float(np.asarray(inside).mean())
gap = 0.90 - coverage
out = ad.AnnData(uns={"metric_ids": ["coverage_at_0.90", "undercoverage_gap"],
                     "metric_values": [coverage, gap]})
out.write_h5ad(par["output"])
