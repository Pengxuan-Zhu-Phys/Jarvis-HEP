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

   The public 2.0.17 build resolves every runtime dependency from public PyPI in a clean environment. It does not require or accept a local Sampler SDK wheel. Use a fresh output directory for each build. Build/test tools are separate from the complete distribution, with no editable install, dependency extras or system-site-packages.

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

## 2.0.17 public release authorization

The user's latest instruction is to publish HEP 2.0.17 while deferring all
Sampler work. This supersedes the earlier candidate-only decision. Build from
public `master` baseline `0bbcf675`, carrying the non-SDK reliability fixes,
CLI agent/YAML assistance, Redis environment commands and their tests/docs.
Both committed SDK integrations `0d859502` / `52876060` and dirty Sampler/Diver
work remain outside this release. The public package has no
`Jarvis-HEP-Sampler` requirement or external `Jarvis sampler` CLI.

Build wheel/sdist from the isolated committed source. Confirm dependency
metadata and package contents exclude the SDK/plugin integration, install all
dependencies from public PyPI, run `pip check`, installed agent/Redis help and
JSON gates, and real Redis/spawn/project-reload smoke. Record fresh source and
installed-package results; the earlier SDK-dependent `1c8c82fc` candidate's
results do not verify this source. Publish through the existing GitHub
Release/Trusted Publishing workflow, then check PyPI and a clean public-index
installation. See [the current release record](docs/release-2.0.17.md).

Current macOS and Linux gates passed: 1426 tests / 477 subtests, 3 skipped,
164 deselected. Linux CI also passed complete installation, pip check, real
runtime/project roundtrip, agent 24 / Redis 39 / release regressions 22.
The [GitHub Release](https://github.com/Pengxuan-Zhu-Phys/Jarvis-HEP/releases/tag/v2.0.17)
and PyPI wheel are public. A fresh complete public-index install passed pip
check, SDK absence, 276 runtime-file identity checks, Redis 39 and agent 24.
Source archive uploads to PyPI did not receive a response; the final bounded request timed out after 36.913 seconds (35-second read timeout). All upload attempts have stopped. The verified source archive is available on GitHub Release; the public wheel and complete public-index installation are verified.

The OIDC publisher job failed with `invalid-publisher`: PyPI had no matching
Trusted Publisher configuration. For this release only, the existing local
PyPI credentials uploaded the exact successful CI artifacts without a rebuild.
This is a documented publication exception; configure the matching Trusted
Publisher before the next normal OIDC release. Public CI artifact hashes differ
from the earlier local validation hashes and are recorded separately in the
[release evidence](docs/release-2.0.17-verification.json).

## Earlier candidate provenance

The earlier SDK-dependent candidate at `1c8c82fc` was verified with a local SDK
and left as a draft. It is superseded for this publication. Its records remain
on the integration branch; their dependency blocker and test counts do not
apply to this master-based source. Current results belong in the
[2.0.17 release record](docs/release-2.0.17.md).

## Sampler requirement

The current source checkout after the SDK-free 2.0.17 package requires
`Jarvis-HEP-Sampler>=1.0.3` without an upper bound. The default environment
cards use the same bound. The unreleased follow-up also requires SDK host
contract 2; public SDK 1.0.3 predates that interface. The maintainer must publish
the matching SDK and update the HEP minimum before releasing this source.
`Jarvis sampler` forwards the SDK `jsampler` commands
and adds HEP's built-in method display. When this HEP source is packaged and
installed, pip resolves a matching sampler release from PyPI. The published
2.0.17 package remains the historical release without this dependency.
