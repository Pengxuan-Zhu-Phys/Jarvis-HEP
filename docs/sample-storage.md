# SAMPLE storage

By default Jarvis keeps a folder of files for every evaluated point under
`outputs/<scan>/SAMPLE/`: calculator input and output files you marked with
`save: true`, and a log for that point. For large scans this can use a lot of
disk space and create millions of small files. If you only need the numbers
stored in `DATABASE/`, turn this off:

```yaml
EnvReqs:
  V2:
    store_samples: false   # default: true
```

The value must be `true` or `false`. You can set it in the project defaults
(`deps/environment_default.yaml`) and override it in a single task card.

## What changes

With `store_samples: false`:

- No `SAMPLE/` folders or `SAMPLE/*.tar.gz` archives are created.
- No per-point logs are written, including the logs of failed points.
- Calculators still run normally. Each point gets its own temporary folder
  (the one `@Sdir` refers to), which is deleted when that point is finished,
  whether it succeeded or failed. Files marked `save: true` are deleted too.
- File paths stored as results in `DATABASE/` point to these temporary files,
  so they will no longer exist after the scan.
- Points that need no files at all (for example, pure Python functions) run
  without any folder.

What is kept:

- All results in `DATABASE/`, including which points failed.
- Resume data, so `Jarvis run TASK.yaml --resume` still works.
- The main logs under `logs/<scan>/`.

This setting takes priority over `sample_directory.enabled`,
`sample_directory.pack`, `archiver.pack_buckets`, `worker.sample_artifacts`,
and `save: true` on individual files.

## Things to know

- The temporary folders are created in your system's temporary directory. On
  Linux and macOS you can choose where by setting the `TMPDIR` environment
  variable, e.g. to a fast local disk on a cluster node.
- If a scan is killed abruptly (for example with `kill -9`, a machine crash,
  or a full disk), some temporary folders may be left behind. Jarvis does not
  clean these up later; delete them by hand if needed.
- Files that a calculator writes to fixed paths outside its point folder are
  not affected.
- Changing this setting never deletes existing `SAMPLE/` folders, and it only
  applies to scans started after the change.
