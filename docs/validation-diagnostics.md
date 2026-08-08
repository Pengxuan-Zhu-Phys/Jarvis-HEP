# Actionable task-card diagnostics

Every task-card validation failure is a stable `JV2-*` error code followed by:

1. the affected YAML path;
2. an explanation of the failed constraint;
3. a `suggestion` that states the next edit; and
4. an `example` when a safe minimal YAML fragment is known.

This also applies before schema validation: malformed task, default-environment,
or runtime-default YAML is reported as `JV2-YAML-001`; a non-mapping task
document is `JV2-YAML-002`; and a missing task file is `JV2-LOAD-001`.
`Jarvis validate --json` exposes the same `suggestion` and `example` fields as
the human-readable output.

`JV2-ENC-001` scans the original parsed task YAML before normalisation and
runtime metadata injection. Keys and string values in that user-authored card
must be ASCII; its path and character positions identify the edit. Generated
values such as `scan_name` and `task_result_dir` are never reported as user
input. Put Chinese explanatory text in a YAML `#` comment instead: comments
are discarded before validation and remain unrestricted. User-provided
non-ASCII diagnostic text is rendered as ASCII `\uXXXX` escapes, including in
the multi-error table.

An encoding failure does not hide independent schema or contract failures:
they are collected in the same report. Only the redundant unknown-key schema
error for the exact non-ASCII key is suppressed, so an ordinary typo beside it
still receives its did-you-mean suggestion.

```text
[error] JV2-SCH-001  $.Calculators.Modules[0].execution.input[0]
        'path' is a required property
        suggestion: Add 'path' at this location using the expected YAML type.
        example:
          name: input
          path: input.json
          type: JSON
```

## Sources of guidance

Schema files may define `x-jarvis-example` next to a card object. The JSON
Schema renderer combines that example with native rule information (`required`,
`enum`, `type`, `additionalProperties`, and numeric/list bounds). This keeps a
format's documentation beside the file that owns its interface.

Python contracts use the same renderer. A path/code-specific guide is used for
sampling, variables, calculator/Opera modules, and V2 runtime configuration.
Every remaining V2 error receives a conservative fallback suggestion, so no
reported validation failure is left without a next action.

Unknown keys in a closed schema object (including the task-card root) also
receive a spelling suggestion when
there is a close allowed field (for example, `naem` suggests `name`). When
multiple unknown keys are reported together, each receives its own spelling
suggestion; very large allowed-key lists are compacted. If YAML
has resolved an unquoted boolean-like scalar where a string is required, the
suggestion explicitly asks for quotes. Accepted numeric strings such as `1e-5`
produce warning `JV2-SCH-003`; they are not rejected for compatibility. An
invalid numeric union reports accepted forms such as `0.05` or `1.0e-5`, rather
than the generic JSON Schema `anyOf` message.

When a card has two or more errors, the human-readable report starts with a
compact `Code` / `YAML path` / `Problem` table, then prints the full actionable
diagnostics. The table explicitly says that it is only a quick reference and
asks the user to consult the Jarvis-HEP YAML settings documentation before
editing the card. It uses ASCII-only borders and escaped text for deterministic
terminal alignment. JSON output remains structured and unchanged.

## Authoring rules

* Examples must be minimal, valid YAML fragments for the failing object.
* Do not embed machine- or project-specific paths in an example.
* Schema format additions must include `x-jarvis-example` and positive and
  negative diagnostic tests.
* Semantic checks that cannot be expressed in JSON Schema must still produce a
  `JV2-*` code through `issue()`; it supplies fallback guidance automatically.
