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


    def manage_directories(self, base_path):
        """
        Manages numbered directories within a given base directory. If the highest-numbered directory contains
        more than 200 files, a new directory with the next number is created.

        Args:
        base_path (str): The path to the base directory to manage.
        """
        # Step 1: Identify all numbered directories
        numbered_dirs = [d for d in os.listdir(base_path) if d.isdigit() and os.path.isdir(os.path.join(base_path, d))]

        if not numbered_dirs:
            # If no numbered directories exist, create the first one and exit
            os.makedirs(os.path.join(base_path, '1'))
            print("Created directory '1' as no numbered directories existed.")
            return os.path.join(base_path, '1')

        # Step 2: Find the highest-numbered directory
        max_number = max([int(num) for num in numbered_dirs])
        max_dir = os.path.join(base_path, str(max_number))

        # Step 3: Check the number of files in the highest-numbered directory
        file_count = len([name for name in os.listdir(max_dir) if os.path.isdir(os.path.join(max_dir, name))])

        # Step 4: Create a new directory if necessary
        if file_count >= 200:
            new_dir_number = max_number + 1
            new_dir_path = os.path.join(base_path, str(new_dir_number))
            os.makedirs(new_dir_path)
            print(f"Created new directory '{new_dir_number}' due to file count greater than 200 in directory '{max_number}'.")
            return new_dir_path
        else:
            print(f"\nDirectory '{max_dir}' contains {file_count} files, no new directory created.\n")
            return max_dir