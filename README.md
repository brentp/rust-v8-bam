# v8bam

This is an experiment with using v8 as a scripting engine. It is more complex than embedding, e.g. lua, but it has great performance (uses v8) and familiar syntax. This has minimal coverage of bam record attributes, but more are easily added.

Rust + V8 bridge to filter BAM/CRAM alignments with JavaScript. A single V8 isolate hosts a reusable `aln` object that is refreshed per record so you can run fast, user-supplied filters.

## General Info

- CLI: `v8bam -e '<js expr>' -o out.bam in.bam` (use `-` for stdin/stdout)
- JS expression can be a boolean expression or a function body; if it lacks `return`, it is wrapped automatically.
- `hasFlag(flag, mask)` is exposed globally for bit tests.
- Multi-threaded BAM I/O via `rust-htslib` thread pool; filter runs single-threaded inside V8.

## JavaScript API (aln object)

- Scalars: `aln.mapq`, `aln.qname`, `aln.flag`, `aln.pos`, `aln.start`, `aln.end`, `aln.chrom`
- Aux tags: `aln.aux("NM")` → number/string/array or `null` if missing
- CIGAR: `aln.cigar` → array of objects `{length, op, consumes_ref, consumes_query}`
  - `op` is one of `Match`, `Ins`, `Del`, `RefSkip`, `SoftClip`, `HardClip`, `Pad`, `Equal`, `Diff`
  - `consumes_ref` and `consumes_query` mirror SAM semantics (e.g., `Match`, `Equal`, `Diff` consume both; `Del`/`RefSkip` only ref; `Ins`/`SoftClip` only query; `HardClip`/`Pad` consume neither)
  - Example: filter out any hard clips
    ```js
    !aln.cigar.some((c) => c.op === "SoftClip");
    ```
- Example filter:
  ```js
  // keep mapped reads with MAPQ>=20, at least one mismatch, and no hard-clipping
  return (
    hasFlag(aln.flag, 0x4) === false &&
    aln.mapq >= 20 &&
    aln.aux("NM") > 0 &&
    !aln.cigar.some((c) => c.op === "HardClip")
  );
  ```

## Using as a Library

- Construct once per expression:
  ```rust
  let mut engine = v8bam::JsBamFilterEngine::new("aln.mapq > 10")?;
  ```
- For each record:

```rust
engine.record_passes(&record, &header_view)?
```

where `record_passes` returns a `Result<bool>`.

- Reuse the same `bam::Record` buffer and header view to minimize allocations.
- The engine owns the V8 isolate/context and reuses a single `aln` object; do not share it across threads without synchronization.
