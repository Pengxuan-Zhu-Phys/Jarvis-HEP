# Releasing Jarvis-HEP

## Package identity

- PyPI distribution: `Jarvis-HEP`
- Python import package: `jarvishep2`
- Command-line entry point: `Jarvis`
- V1 releases stop at the `1.x` series; V2 starts at `2.0.0` on the same
  `Jarvis-HEP` PyPI project.

PyPI normalizes distribution names, but it does not infer them from the import
directory. Keep `[project].name = "Jarvis-HEP"` in `pyproject.toml` even though
the source directory is named `jarvishep2`.

## One-time Trusted Publishing setup

The release workflow uses PyPI Trusted Publishing (GitHub OIDC), so no PyPI API
token is stored in GitHub.

1. In GitHub, open **Settings → Environments**, create an environment named
   `pypi`, and require a maintainer's approval for deployment if the repository
   plan supports it.
2. In the PyPI `Jarvis-HEP` project, open **Manage → Publishing**, add a GitHub
   Actions trusted publisher, and enter:

   | Field | Value |
   | --- | --- |
   | Owner | `Pengxuan-Zhu-Phys` |
   | Repository | `Jarvis-HEP` |
   | Workflow filename | `publish-pypi.yml` |
   | Environment | `pypi` |

The workflow file must be merged into the repository's default branch before
the first release is published.

## Release procedure

1. Complete the development ledger and run the fast gate plus relevant runtime tests. Add intended new files to Git; the helper snapshots tracked files and their current edits. Keep unfinished packages out of the release snapshot.
2. Run `./jhrel <version> --output dist/<version>` to update both source version fields, build wheel and sdist, install the complete wheel with all dependencies in a clean venv, run `pip check`, installed CLI/runtime and agent-authoring smoke, and the 2.0.15 regression contracts. `--dry` is an alias for this default; neither uploads. The helper emits `release-verification.json` and `SHA256SUMS` only after verification. `--install` is not supported.

   For unpublished dependency candidates use `--sampler-dist DIR` and `--operas-dist DIR`. These install the dependency wheels with every declared dependency; publish dependencies before releasing HEP. Use a fresh output directory for each build. Build/test tools are separate from the complete distribution, with no editable install or dependency extras.

3. Merge the release commit into `master`.
4. Create a GitHub Release whose tag is exactly `v<version>`, for example:

   ```bash
   gh release create v2.0.0 --target master --generate-notes --title "Jarvis-HEP 2.0.0"
   ```

Publishing the GitHub Release triggers `.github/workflows/publish-pypi.yml`.
The workflow checks that the tag, both source version declarations, V2 major
version, distribution name, and wheel contents agree. It then runs the default
test suite, builds one wheel and one source distribution, and publishes those
exact artifacts to PyPI using a short-lived OIDC credential.

Do not upload with a local token as part of the normal release process. Do not
enable `skip-existing`: a duplicate version should fail visibly instead of
hiding a release mistake.

## Upgrading from V1

Existing V1 users can upgrade directly from `Jarvis-HEP` 1.x to 2.x, but V2 is
a major-version migration: the old `jarvishep` import package is replaced by
`jarvishep2`.

## 2.0.15 handoff

See [the development ledger](docs/development-ledger-2.0.15.md) and [release notes](docs/bugfix-2.0.15/release-notes.md). The agent prepares verified candidates locally; the user creates the final release/tag and publishes. No upload, tag creation or push is performed during candidate preparation.

## 2.0.16 publication authorization

The user authorized committing and publishing the CLI agent improvements on 2026-10-13. This supersedes the earlier candidate-only handoff for this release. See [the 2.0.16 release record](docs/release-2.0.16.md). Build from a clean release commit checkout to exclude unrelated working-tree changes, and publish the required Sampler SDK before HEP.

## 2.0.17 Redis publication authorization

On 2026-10-15 the user authorized committing the Redis improvements and
publishing a new HEP version, while explicitly excluding the current Sampler
changes. See [the 2.0.17 release record](docs/release-2.0.17.md). Include the
Redis environment/command-group implementation, its tests and documentation,
and the two HEP version declarations. Preserve unrelated Sampler/runtime,
Diver, package-data and review-file edits outside the commit and candidate.

Build wheel and sdist from the clean committed source, not the dirty working
tree. Verify their contents, clean complete installation, `pip check`, the
installed Redis parent/child help and JSON/agent behavior, and real Redis/spawn
runtime smoke. Earlier Redis test results are historical evidence until this
release snapshot passes its own checks.

The existing committed baseline still requires `Jarvis-HEP-Sampler>=0.1.0,<0.2`.
An already verified local SDK 0.1.0 wheel may be used unchanged for candidate
validation; that does not prove installation from public indexes. Public HEP
release remains blocked while this required package is absent from PyPI.
Prepare a draft and record the blocker; do not modify or publish the SDK under
this Redis-only authorization. Follow the normal GitHub Release/Trusted
Publishing path when the dependency requirement is resolved.

Later in the same session, the user selected candidate preparation only and
retained SDK publication as a follow-up. Complete the Redis-only commits,
candidate verification and draft delivery; public HEP/SDK publication remains
follow-up work.
