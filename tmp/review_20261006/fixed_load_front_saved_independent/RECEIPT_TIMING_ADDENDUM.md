# Receipt timing addendum

When this independent cached-array review was completed, `run01/result.json` and the experiment README were available, but the run folder did not yet contain `outer_process_receipt.json` or `OUTER_PROCESS_ADDENDUM.md`. The original `execution_receipt.json` therefore treated the author-reported process status as README-reported evidence and recorded that raw wrapper streams were unavailable.

Those two files were added afterward. The structured receipt is 4,603 bytes, SHA-256 `2A0615A0CBFAEE065669BA28CF56E0D979811ACBFD377E790F65FBFD3F3E8386`; the Markdown addendum is 2,929 bytes, SHA-256 `D12F5752FEFCC702016E9A215FE8A6B9FBF9F159BFE15CC35FE681559A85AC9D`. They record a post-hoc reconstruction from the existing tool response and saved result: wrapper-reported `rc=0`, elapsed `0.5507252 s`, no timeout, and result SHA-256 `a60e51c85c26510c358ee57cafa4bfdda6ed89ea463d99bcc86551f64220732c`. Both explicitly state that raw stdout/stderr bytes were not persisted and that no execution was repeated to create the receipt.

This later metadata improves the timing/provenance record; it is not a contemporaneous raw process capture and does not change the cached-array audit, its interpretation, or the recorded result. No numerical work was repeated for this addendum.
