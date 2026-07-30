# Actionable task-card diagnostics

Every task-card validation failure is a stable `JV2-*` error code followed by:

1. the affected YAML path;
2. an explanation of the failed constraint;
3. a `suggestion` that states the next edit; and
4. an `example` when a safe minimal YAML fragment is known.

This also applies before schema validation: malformed YAML syntax is reported as
`JV2-YAML-001`, a non-mapping document as `JV2-YAML-002`, and a missing task
file as `JV2-LOAD-001`. `Jarvis2 validate --json` exposes the same `suggestion`
and `example` fields as the human-readable output.

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

## Authoring rules

* Examples must be minimal, valid YAML fragments for the failing object.
* Do not embed machine- or project-specific paths in an example.
* Schema format additions must include `x-jarvis-example` and positive and
  negative diagnostic tests.
* Semantic checks that cannot be expressed in JSON Schema must still produce a
  `JV2-*` code through `issue()`; it supplies fallback guidance automatically.
