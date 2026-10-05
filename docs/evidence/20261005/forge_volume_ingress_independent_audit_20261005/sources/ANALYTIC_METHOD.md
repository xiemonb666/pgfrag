# Independent analytic reference method — frozen annular-cylinder ROI

This is the proposed independent numerical oracle and data-adapter plan for the actual-output audit. The production worker, its CAD booleans, its reference ROI summary, and its volume report are not used to compute expected geometry. The expected values come only from the fixture geometry stated by root: a straight annulus with outer radius 10 mm, inner radius 8 mm, z=[0,10] mm, and the fixed ROI x=[7.5,10.5], y=[-1,1], z=[4,6] mm.

## Analytic construction

For each XY grid rectangle, compute the area and first moments of a disk of radius R clipped to the rectangle, then subtract the corresponding radius-8 disk from the radius-10 disk. This is valid because the annular cross-section is exactly the outer disk minus its concentric inner disk.

At a fixed x, the circle arc is s(x)=sqrt(R^2-x^2). Clip the vertical interval [−s,s] against each rectangle [y0,y1]. Split the x integration interval at every x where s(x)=|y0| or |y1|; between these breakpoints the active upper and lower bounds are fixed as either a constant rectangle edge or ±s(x). Integrate each segment with closed-form primitives:

- ∫s dx = 0.5*(x*s + R^2*asin(x/R))
- ∫x*s dx = −(R^2−x^2)^(3/2)/3
- ∫0.5*s^2 dx = 0.5*(R^2*x − x^3/3)

For the active upper/lower boundaries yU and yL, the segment contributions are area ∫(yU−yL)dx, x first moment ∫x(yU−yL)dx, and y first moment 0.5∫(yU^2−yL^2)dx. The resulting XY centroid is first-moment/area. Because each z slab is uniform and fully inside the cylinder height, its volume is area×Δz and z centroid is the slab midpoint.

The script analytic_roi_reference.py implements only these standard-library formulas; it imports no CAD package and does not execute on import. It generates the deterministic 96 rows in k,j,i order and computes the whole ROI directly, separately from the sum of the 96 cell integrals. The independent whole-ROI COM is also compared with the volume-weighted COM from those cell expectations.

## Planned read-only adapter after output freeze

After the writer publishes and root pins the actual-output freeze, adapt only to the frozen artifact names and real CSV columns. Verify each used artifact SHA against that freeze before parsing. Read raw rows, not the producer summary. Cross-check the CSV against the cell ledger by explicit cell index/ijk and component/point IDs; if those mappings or source classifier evidence are absent, report them as UNKNOWN rather than reconstructing them from summary claims.

The adapter will independently aggregate delivered point volumes and volume-weighted centroids per grid cell, compare those against the analytic cells, and separately compare the analytic direct-ROI integral against its 96-cell sum. It will verify any supplied mm³-to-m³ conversion by multiplying by 1e−9, and mass by rho×volume in SI using only the pinned declared rho=2100 kg/m³. It will report absolute and relative deltas without promoting agreement to calibration or physical validation. Void-only cells must have analytic zero volume and no delivered material point; positive cells must map to the frozen cell bounds. Centroid classifier evidence must be present in the raw ledger and be IN for every emitted point; missing evidence remains UNKNOWN.

The two planned images use only verified actual rows and the independent analytic table:
1. Equal-aspect XY cell map with analytic material fill per XY cell, actual volume-centroid markers, ROI edges, and inner/outer circles.
2. Per-cell analytic versus delivered volume, declared-density mass check, and signed volume delta; cells are labeled by stable grid index/ijk and zero cells remain visible.

Root should review the adapter once its exact input schema is known. After that review only, run one read-only audit bounded to 60 seconds. Do not invoke CadQuery, OpenCascade, the production worker, or a solver from the audit.