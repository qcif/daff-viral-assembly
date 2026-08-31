# Refdata blob download speed test — progress

Goal: measure azcopy download throughput from **standard** vs **premium**
blob storage to an Azure Batch node, to decide which tier to use for
staging the ~830GB wf4 reference data.

## Fixtures

- Test containers: `test-refdata` on both `daffstandard` and `daffpremium`
  (throwaway — delete when done).
- Test blob (planned): `k2_core_nt_20251015.tar.gz` — uploaded from the
  remote VM where the refdata lives.
- Speedtest pool: `view_speedtest` — 1× `Standard_L8as_v3` dedicated node
  (smallest allowed L-series SKU on this Batch account; `L4as_v3` was
  rejected as too small). One task slot per node.

## Files in this directory

| File | Purpose |
|---|---|
| `refdata-speedtest-upload.sh` | Run on the remote VM: uploads a file to `test-refdata` in both storage accounts via azcopy (generates short-lived write SAS from `.env.azure` account keys). |
| `pool-setup-speedtest.json.template` | Pool config template; `<SETUP_SCRIPT_SAS_TOKEN>` filled at render time. |
| `setup-speedtest.sh.template` | Start-task template. Placeholders: `<TEST_BLOB_NAME>`, `<STD_SAS_TOKEN>`, `<PREM_SAS_TOKEN>`. Comment/uncomment the `SOURCE_URL` line to pick which tier to benchmark. |
| `rendered/setup-wf4-speedtest.sh` | Rendered start task, uploaded to `daffstandard/scripts/setup-wf4-speedtest.sh`. |
| `rendered/pool-setup-speedtest.json` | Rendered pool config used to create the pool. |

## Status

- [x] Created `test-refdata` container on `daffstandard` (2026-07-21).
- [x] Created `test-refdata` container on `daffpremium` (2026-07-21).
- [x] Generated 30-day read+list SAS tokens for both `test-refdata` containers.
- [x] Rendered `setup-wf4-speedtest.sh` — currently benchmarks **STANDARD** tier;
      premium `SOURCE_URL` line is commented out.
- [x] Uploaded rendered start-task script to `daffstandard/scripts/setup-wf4-speedtest.sh`
      and generated a 30-day read SAS.
- [x] Created Batch pool `view_speedtest` (1× `Standard_L8as_v3`, 1 dedicated node).
- [ ] Upload `k2_core_nt_20251015.tar.gz` to both `test-refdata` containers
      from the remote VM using `refdata-speedtest-upload.sh`.
- [ ] Wait for pool node to come up and start task to run; capture standard-tier
      throughput from `az_node_logs view_speedtest stderr`.
- [ ] Re-render `setup-wf4-speedtest.sh` with the premium `SOURCE_URL` line
      un-commented (and standard line commented). Re-upload to scripts container.
- [ ] Recreate node so new start task runs: `az_pool_resize 0 view_speedtest --yes`
      then `az_pool_resize 1 view_speedtest --yes`.
- [ ] Capture premium-tier throughput from start-task logs.
- [ ] Record results table below and decide on tier.
- [ ] Teardown: `az_pool_delete view_speedtest`; delete `test-refdata` containers
      on both accounts; delete `scripts/setup-wf4-speedtest.sh`.

## Results

| Tier | File | Size | Elapsed (s) | Speed (Mbps) | Node SKU | Date |
|---|---|---|---|---|---|---|
| standard | | | | | Standard_L8as_v3 | |
| premium  | | | | | Standard_L8as_v3 | |

## Notes

- SAS tokens expire 2026-08-20 (30 days from 2026-07-21). Regenerate if
  the test slips past that.
- `L4as_v3` (originally requested) was rejected by the Batch account —
  minimum allowed L-series SKU is `L8as_v3`. NIC bandwidth is higher on
  `L8as_v3` (~12.5 Gbps) than `L4as_v3` (~6.25 Gbps) but still well below
  prod `L48as_v3` (~40 Gbps), so absolute throughput will still under-read
  vs prod — the standard-vs-premium *ratio* is the useful signal here.
