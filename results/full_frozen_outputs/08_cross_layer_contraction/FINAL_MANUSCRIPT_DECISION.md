# Final manuscript decision after cross-layer discrimination analysis

Decision: PIVOT — continue only as a structural-calibration / evidence-boundary methodological paper; do not continue mechanism-rescue analyses.

Observed cross-layer GJSD:
- Herb: 0.5514179554
- Target: 0.0527603501 (9.568% of herb)
- GO-BP: 0.0042186265 (7.996% of target; 0.765% of herb)
- KEGG: 0.0060162013 (11.403% of target; 1.091% of herb)

Exact Curveball structural null (10,000; exact herb degree, target degree, and edge count preservation):
- Target: null mean 0.04014653, Z=3.77965, upper P=0.0028997.
- GO-BP: null mean 0.00319976, Z=1.31309, upper P=0.0894911.
- KEGG: null mean 0.00324114, Z=2.52331, upper P=0.0232977.

Retention-ratio null:
- Target/herb observed 0.09568; upper P=0.00290.
- GO/target observed 0.07996; upper P=0.42116.
- KEGG/target observed 0.11403; upper P=0.09719.

Interpretation:
1. There is a large, monotonic projection-associated contraction of cross-syndrome discrimination from herbs to targets to functions.
2. Real herb-target assignments retain more target-space differentiation than exact degree structure alone.
3. Functional differentiation is weak/diffuse. GO does not exceed the absolute structural null. KEGG shows modest global differentiation above its absolute structural null, but target->KEGG retention is not unusually high relative to the exact structural-null retention distribution and the frozen pathway-level analysis found 0/621 FDR-supported syndrome-specific pathway tests.
4. Therefore the paper must not claim that topology fully explains functional specificity loss, nor that syndrome-specific pathways have been established.
5. Combined with CAP, psoriasis, URTI, and DCSMI results, the defensible central message is an evidence-calibration boundary: treatment-derived representations can remain computationally differentiated in target space while progressively losing interpretable specificity and failing independent patient-level validation.

Investment recommendation:
- Continue to one manuscript rewrite cycle centered on structural calibration, projection-associated contraction, and prospective falsification.
- Do not add new enrichment databases, diseases, machine-learning models, or alternative cutoffs to seek positive mechanisms.
- If peer-review fit is poor after targeted journal selection / pre-submission assessment, pivot the project rather than initiating another analysis-rescue cycle.
