# Provenance scope addendum

The v2 audit result is `audit_result_v2.json` (SHA-256 `930495dc9b4f8ed1c7a0a0bd5b496712cc883c4a3477765ad93137f1e8bfc8e6`). Its external execution record is `execution_record_v2.json` (SHA-256 `818eeb060c3136a92a8d35df222db5d3df682afcd47205674f9e4770ab751312`).

The embedded v1 field `successful_v1_launch.v1_receipt_sha256` contains the SHA-256 of the original operator `run_record.json` (`101f948ea40c40b63845c60df417efad84343d9fd41ad2ebf875a83cbab85f2a`). It is not a hash of a separate external v1 launch receipt. The v1 launch commands, exit codes, and wall times in that record are provenance copied from the tool execution history.

This note clarifies provenance only. The v1 and v2 source and result artifacts are unchanged.
