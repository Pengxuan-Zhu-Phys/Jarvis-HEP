import subprocess

class Package:
    def __init__(self, name, required_modules, installed, installation):
        self.name = name
        self.required_modules = required_modules
        self.installed = installed
        self.installation = installation

    def install(self, library_settings):
        if self.installed:
            print(f"{self.name} is already installed.")
            return
        for command in self.installation["commands"]:
            formatted_command = command.format(**library_settings, **self.installation)
            subprocess.run(formatted_command, shell=True, check=True)
        self.installed = True
        print(f"{self.name} installation complete.")

class Library:
    def __init__(self, settings):
        self.settings = settings
        self.packages = []

    def add_package(self, package):
        self.packages.append(package)

    def install_packages(self):
        for package in self.packages:
            if not package.installed:
                if package.required_modules:
                    for dependency_name in package.required_modules:
                        dependency = next((p for p in self.packages if p.name == dependency_name), None)
                        if dependency and not dependency.installed:
                            dependency.install(self.settings)
                package.install(self.settings)

# Example usage
settings = {
    "package_path": "/External/Library",
    "make_paraller": "16",
    "root_path": "/path/to/root"
}

library = Library(settings)
library.add_package(Package("Delphes", [], False, {"path": "/External/Library/Delphes", "source": "/External/Library/Source/Delphes-3.5.0.tar.gz", "commands": [...] }))
library.add_package(Package("HepMC", [], False, {"path": "/External/Library/HepMC", "source": "/External/Library/Source/HepMC-2.06.09.tar.gz", "commands": [...] }))
library.add_package(Package("Pythia8", ["HepMC"], False, {"path": "/External/Library/Pythia8", "source": "/External/Library/Source/pythia8230.tgz", "commands": [...] }))

library.install_packages()
