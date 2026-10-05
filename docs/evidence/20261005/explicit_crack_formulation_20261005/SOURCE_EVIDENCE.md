# Source evidence and extraction record

Scope: research and design for the selected explicit 3D CPM/EFGM LEFM branch. No solver was run and no tests or production files were changed.

## Primary files

1. **Ai & Augarde (2018), accepted article.** “An adaptive cracking particle method providing explicit and accurate description of 3D crack surfaces,” *International Journal for Numerical Methods in Engineering* 114(12), 1291–1309, DOI [10.1002/nme.5786](https://doi.org/10.1002/nme.5786). The official record is [Durham Research Online](https://durham-repository.worktribe.com/output/1335971/an-adaptive-cracking-particle-method-providing-explicit-and-accurate-description-of-3d-crack-surfaces); the supplied file endpoint `https://durham-repository.worktribe.com/OutputFile/1335985` returned a Cloudflare 403 challenge in this environment. The browser PDF parser could read a 21-page accepted-version copy at `https://scispace.com/pdf/an-adaptive-cracking-particle-method-providing-explicit-and-2ki0z9yv6s.pdf`, whose first page identifies it as a Durham Research Online accepted version. A direct shell request to that mirror returned HTTP 202 `text/html` with zero bytes. **Accepted-article local SHA-256: UNKNOWN (source bytes were not obtainable).** Do not infer a hash from rendered/text-extracted browser content.

   Useful printed journal locators: §2 / Eqs. (1)–(6), pp. 1293–1295; §3.1 surface representation and visibility, p. 1296; §3.2 SIF interaction integral, pp. 1297–1298; §3.3 growth rule, pp. 1298–1299; Algorithm 1, p. 1300; §4.2 inclined penny crack and analytical SIFs, pp. 1301–1302, Eqs. (20)–(22). These locators refer to the article’s printed pagination; browser PDF page indices differ.

2. **Ai (2018), author’s official doctoral thesis.** *Computational fracture modelling by an adaptive cracking particle method*, Durham University, [official PDF](https://etheses.dur.ac.uk/12679/1/WAi_PhDthesis.pdf). Downloaded as `ai_phd_thesis_official.pdf`, 20,667,044 bytes; SHA-256:

   ```text
   51F7480AAECC8E6BA68B148447BE595D5A68E50FB3D33A2A184B377B586740F0
   ```

   Extracted page text is preserved in `thesis_ch2_selected_pages.txt` and `thesis_ch6_selected_pages.txt`. The thesis is the author’s primary supporting source for missing article prose and benchmark parameters, but article-specific differences are recorded rather than blended away.

## Formula/page map

| Topic | Source locator | Use and limits |
|---|---|---|
| MLS basis, support, kernel example, shape functions and gradients | Thesis Ch. 2 §2.2, printed pp. 32–33; PDF pp. 57–58; Eqs. (2.1)–(2.10) | Explicit 3D linear basis `[1,x,y,z]ᵀ`; MLS moment matrix and analytic gradient. Compact quartic spline is presented as an example; support radius is said to be typically 2–4 average nodal spacings. Neither exact kernel nor exact support factor for 2018 §4.2 is stated. |
| Standalone 3D MLS primitive | [Experiment README](../../../experiments/mls3d_kernel_20261005/README.md), [primitive contract](../mls3d_kernel_20261005/CONTRACT.md), [source review](../mls3d_kernel_20261005/SOURCE_REVIEW.md) | The root-selected primitive uses Wendland C2 with explicit support configuration as an engineering choice. The author's quartic spline remains a separate source example; neither setting is claimed for the 2018 §4.2 benchmark. The primitive is `MATHEMATICAL_PRIMITIVE_ONLY`, not a coupled crack backend. |
| Surface visibility | Article §3.1, p. 1296; thesis Ch. 6 §6.2, printed pp. 128–130; PDF pp. 153–155, Fig. 6.1–6.3 | Finite crack patch intersection controls visibility; triangular facets define nonplanar local support segments. Crack-front particles are not split. Article extraction says copy offset `0.01 × support`; thesis says `0.001 × support` as an example. |
| Crack opening | Article Eq. (2), p. 1294 | Difference of the MLS traces interpolated from `S+` and `S−` side supports. |
| SIF integration and local frame | Article §3.2, pp. 1297–1298, Eqs. (7)–(9); thesis Ch. 6 §6.3.1, printed pp. 131–132; PDF pp. 156–157, Eqs. (6.1)–(6.3), Fig. 6.4–6.5 | Box `2c × 2c × 2lc`, product `q` weight and auxiliary fields. Numeric `c,lc` are not reported. Source Eq. (8)/Eq. (6.2) prints `+ μ KIII KIIIaux`, conflicting with its own energy-consistent Eq. (12)/(6.4) scaling; preserve as `SOURCE_FORMULA_CONFLICT`. Independent backend derivation uses `1/μ` and requires a pure-mode-III analytic check; it is not attributed to the paper. |
| Direction, equivalent SIF, front increment | Article §3.3, pp. 1298–1299, Eqs. (10)–(14); thesis Ch. 6 §§6.3.2–6.3.3, printed pp. 132–134; PDF pp. 157–159, Eqs. (6.4)–(6.6), Algorithm 6.1 | Maximum circumferential stress angle ignores mode III for direction; Paris-law-scaled front increment; `m=3.32` is specified in the examples. Numeric `Δa_max` is not stated. Quad strips are divided into three triangles via a midpoint; the midpoint is not advanced in the simple algorithm. |
| General 3D example discretization and rigid constraints | Thesis Ch. 6 §6.4, printed pp. 135–136; PDF pp. 160–161, Eqs. (6.7)–(6.9) | `2w=200 mm`, initial regular `15³` points, volume `14³` cells × `4³` Gauss, loaded surface `14²` cells × `8²` Gauss. Six Lagrange multiplier constraints set mean translation and rotation to zero. |
| Inclined penny crack | Article §4.2, pp. 1301–1302, Eqs. (20)–(22); thesis Ch. 6 §6.4.2, printed pp. 139–141; PDF pp. 164–166, Fig. 6.12–6.15 and Eqs. (6.11a–c) | `a=w/5`, `β=π/6` from x-axis within xz plane, same baseline `E=100 GPa`, `ν=.3`, `σ=100 MPa`, target error `ηt=.02`. The thesis explicitly supplies load units; accepted-article text extraction renders `N/m²`, recorded as a source-text discrepancy pending visual article verification. Analytic SIFs are infinite-body values; source reports finite-cube end-angle deviations. |

## Bounded extraction and decisions

- The official thesis PDF is the only source file downloaded to this directory. Its full file hash is recorded above. A text-only extraction of the author’s accepted paper was available through the browser tool, but no accepted-paper bytes were retained and no hash is asserted.
- The selected §4.2 benchmark does **not** specify exact MLS kernel/support factor, interaction-box `c,lc`, `Δa_max`, or a numeric tolerance for matching the SIF curves. These remain `UNKNOWN` for an exact paper-reproduction claim.
- A dimensionless static linear-elastic fixture with `E=1`, `σ=1e-3`, and `ν=.3` preserves the paper-supported ratio `σ/E=1e-3` while testing mathematical mechanics primitives. It is not a paper reproduction, calibrated material fracture, or GT. The current research contract and MLS primitive are `NOT_COUPLED`; there is no mechanics-grounded surface or model R2. A future result may qualify as `SIMULATED_GEOMETRIC_GT` only after explicit surface/mechanics coupling and surface-consistency plus material-partition gates pass; `MATERIAL_UNCALIBRATED` remains until material validation.
- The 2019 Ai et al. 2D CPM article (DOI [10.1016/j.enganabound.2019.03.008](https://doi.org/10.1016/j.enganabound.2019.03.008)) gives `d_s=2.2` for its own 2D examples. This is not evidence that the 3D §4.2 benchmark used that factor.
