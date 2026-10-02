Target Projection v1.1 — FROZEN
================================

Status
------
FROZEN_UPSTREAM_BEFORE_CAP_ALIGNMENT

Primary downstream seed
-----------------------
LOPO_FROZEN_primary_target_weights_for_RWR_v1_1.csv
Use column: ubiquity_corrected_sum1
Group by: syndrome_code
Target identifier: target_id

Frozen target-layer processed-form mappings
-------------------------------------------
制陈皮 -> 陈皮
炙百部 -> 百部
焦山楂 -> 山楂

Prescription-layer herb identities are preserved. These parent mappings are
used only to obtain molecular targets.

Frozen projection
-----------------
For syndrome s and target t:
  degree-normalized contribution = W_hs / |T_h|
  IDF_t = log((N_H + 1)/(n_t + 1)) + 1
  target score = sum_h [W_hs / |T_h| * IDF_t]
  final vector = target score / sum_t(target score)

N_H is the number of unique herbs with at least one valid target in the
completed study target table. Processed-form aliases do not expand this IDF
reference universe.

Validated frozen run
--------------------
Target DB herbs with >=1 valid target: 172
Unique targets: 374
Unique herb-target pairs: 13,285
Degree-normalization max mass error: 1.1102230246251565e-16

LOPO mapping coverage:
  PDOL 85/95 herbs; weight fraction 0.940547
  PHOL 107/121 herbs; weight fraction 0.900780
  WHIL 87/96 herbs; weight fraction 0.893701

LOPO ubiquity-corrected Spearman similarity:
  PDOL-PHOL 0.928626
  PDOL-WHIL 0.881911
  PHOL-WHIL 0.936850

No CAP transcriptomic data were used to select, tune, or freeze this projection.
Any later change to mapping rules, target database, IDF universe, or projection
formula requires a new version; v1.1 must not be overwritten.
