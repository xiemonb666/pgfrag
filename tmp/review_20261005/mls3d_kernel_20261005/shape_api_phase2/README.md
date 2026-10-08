# MLS shape API phase 2

This review bundle adds a bounded direct MLS shape-function interface while
preserving the already reviewed and tested phase1 sources. See [CONTRACT.md](CONTRACT.md)
for the exact API/schema and [SOURCE_REVIEW.md](SOURCE_REVIEW.md) for the scope
and review gates. All seven focused phase2 acceptance tests passed; receipts
are in `phase2_focused_test/`.

All claims remain `MATHEMATICAL_PRIMITIVE_ONLY` / `NOT_COUPLED`; this does not
create fracture geometry or R2 supervision.
