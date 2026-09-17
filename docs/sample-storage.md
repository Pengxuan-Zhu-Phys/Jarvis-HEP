# SAMPLE storage policy

Implemented 2026-09-17: `EnvReqs.V2.store_samples` is a strict YAML boolean,
defaulting to `true`. Project environment defaults can be overridden by a task.

```yaml
EnvReqs:
  V2:
    store_samples: false
```

When false, the runtime does not create SAMPLE or allocate numbered buckets,
does not tar buckets, and does not open or buffer per-sample logs (including
failure details and calculator installation reuse logs). Process-level logs,
errors, DATABASE observables/status and durable resume indices remain enabled.
The setting overrides `sample_directory.enabled`, `sample_directory.pack`,
`archiver.pack_buckets`, persistent IO `save: true` retention, and `worker.sample_artifacts`.

Pure numeric samples need no directory. Calculators and workflows using `@Sdir`
receive an isolated directory from Python's system temporary-directory policy
(`TMPDIR` on Unix can select local scratch). IO files, including copies needed
after a calculator PackID is released, remain available to downstream steps
until that sample completes. The Worker cleans its owned scratch directory in
the finalization path, also on computation or archive-submission failure. No
scratch directory is sent to the Archiver as a persistent product directory.
File-path observables are temporary references, not retained downloadable files.

Existing SAMPLE directories are never deleted by changing this setting. YAML
changes apply to newly started runs; they do not reconfigure running processes.
Required calculator files still incur temporary IO, and calculator-authored
files at explicit paths outside SAMPLE remain under that calculator's control.
SIGKILL, machine failure or filesystem errors can leave scratch behind; this
setting does not add a background scratch scavenger. Installation metadata and
process logs are retained for operational diagnostics.

## Implementation and verification

The normalized runtime flag feeds Worker blueprints, bucket configuration and
Archiver policy. Sample owns scratch allocation, log suppression and cleanup;
the existing IO copy policy runs inside scratch to preserve downstream readers.
On-demand `@Sdir` resolution updates the shared sample info so repeated resolution
uses the same directory and cleanup knows which directory it owns. DATABASE
batch commit and resume behavior are unchanged.

Regression coverage: strict boolean validation, inherited defaults and task
overrides, default-on behavior, lazy `@Sdir` reuse, silent bound loggers, numeric
samples without directories, scratch cleanup on success/failure, DATABASE
status and resume prefix, and a real calculator check against golden results.

This change uses the current implementation and Jarvis-Books V2 design docs as
authorized in the implementation conversation: the five historical canonical
documents named by the workspace AGENTS.md are absent from this checkout.
