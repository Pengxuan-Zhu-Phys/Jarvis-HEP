"""Safe YAML loading with explicit duplicate-key diagnostics.

Check original nodes before SafeLoader flattens merge keys: inherited values
may legitimately be overridden, but two explicit keys must never disappear.
"""

from __future__ import annotations

from typing import Any

import yaml

_MERGE_KEY = object()


class DuplicateYamlKeyError(yaml.YAMLError):
    def __init__(self, key: Any, path: str, first_mark: Any, second_mark: Any):
        self.key = key
        self.path = path
        self.first_mark = first_mark
        self.problem_mark = second_mark
        super().__init__(f"duplicate YAML key {key!r} at {path}")


class UniqueKeySafeLoader(yaml.SafeLoader):
    def get_single_data(self) -> Any:
        node = self.get_single_node()
        if node is None:
            return None
        self._check_keys(node, "$", set())
        return self.construct_document(node)

    def _check_keys(self, node: Any, path: str, visited: set[int]) -> None:
        # Aliases may refer back to the same node. Visit each original mapping
        # once; this also keeps recursive aliases from looping in this check.
        if id(node) in visited:
            return
        visited.add(id(node))
        if isinstance(node, yaml.MappingNode):
            seen: dict[Any, Any] = {}
            for key_node, value_node in node.value:
                is_merge = key_node.tag == "tag:yaml.org,2002:merge"
                key = "<<" if is_merge else self.construct_object(key_node, deep=True)
                identity = _MERGE_KEY if is_merge else key
                try:
                    first = seen.get(identity)
                    if first is not None:
                        raise DuplicateYamlKeyError(key, f"{path}.{key}", first, key_node.start_mark)
                    seen[identity] = key_node.start_mark
                except TypeError as exc:
                    raise yaml.constructor.ConstructorError(
                        "while constructing a mapping", node.start_mark,
                        "found unhashable key", key_node.start_mark,
                    ) from exc
                self._check_keys(value_node, f"{path}.{key}", visited)
        elif isinstance(node, yaml.SequenceNode):
            for index, child in enumerate(node.value):
                self._check_keys(child, f"{path}[{index}]", visited)
