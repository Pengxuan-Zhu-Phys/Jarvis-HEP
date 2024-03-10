#!/usr/bin/env python3
import os, sys 
import json

class Base():
    def __init__(self) -> None:
        self.path: Dict[str, str]   = {
            "jpath":    os.path.abspath(os.path.join(os.path.dirname(__file__), "..")) 
        }
        self.path['args_info'] = self.decode_path("&J/src/card/argparser.json")
        self.path['logger_config_path'] = self.decode_path("&J/src/card/jarvis_logging_config.yaml")
        self.load_schema_paths()

    def decode_path(self, path) -> None: 
        """
        Resolves special markers in the provided path.

        Parameters:
        - path: The path to be resolved.
        - jarvis_root: The root directory path of Jarvis.

        Returns:
        - The resolved full path.
        """
        # Replace the Jarvis root directory marker &J
        if "&J" in path:
            path = path.replace("&J", self.path['jpath'])

        # Replace the user home directory marker ~
        if "~" in path:
            path = os.path.expanduser(path)
        
        return path  

    def load_schema_paths(self) -> None:
        self.path['preference'] = self.decode_path("&J/src/card/preference.json")
        with open(self.path['preference'], 'r') as f1:
            js = json.loads(f1.read())
            for kk, vv in js["Schema"].items():
                js["Schema"][kk] = self.decode_path(vv)
            
            self.path.update(js['Schema'])
            self.path['logo'] = self.decode_path(js['logo'])
