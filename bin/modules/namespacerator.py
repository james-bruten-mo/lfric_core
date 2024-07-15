###############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
###############################################################################
"""
Implements namespacing tool.
"""
from abc import ABC, abstractmethod
import re
from pathlib import Path
from re import compile as re_compile
from subprocess import run
from typing import Callable, Dict, List, Optional

IFORT_LIMIT = 63


def error_callback(symbol: str) -> str:
    """
    Raises an exception for over long symbol names.

    Does not attempt to remidy the situation.
    """
    raise Exception(f"Symbol '{symbol}' too long.")


class AddNamespace:
    """
    Mutates Fortran source by adding a namespace string to module names.
    """
    def __init__(self, name: str,
                 long_symbol_callback: Callable[[str], str]
                 = error_callback):
        """
        :param name: Namespace name.
        :param long_symbol_collback: Called if the generated symbol is too
                                     long. Defaults to raising an exception.
        """
        self.__long_symbol_callback = long_symbol_callback

        # Ensure the namespace name ends with an underscore.
        #
        if name.endswith('_'):
            self.__name = name
        else:
            self.__name = name + '_'

        # Keep track of renames
        #
        self.__renames: Dict[str, str] = {}

        # Compile RegExes for later use.
        #
        self.__module_pattern = re_compile(
            r'(?P<head>\s*(end\s*)?module\s*(?!procedure))'
            r'(?P<name>\w+)(?P<tail>.*)$', re.IGNORECASE
        )

    def rename_module(self, source: str) -> str:
        """
        Adds namespace to Fortran module source.

        :param source: Fortran source code.
        :return: Modified Fortran source code.
        """
        new_text: List[str] = []
        current_module: Optional[str] = None
        for line in source.splitlines(keepends=True):
            match = self.__module_pattern.match(line)
            if match:
                current_module = match.group('name')
                name = current_module
                for _ in range(5):
                    new_name = self.__name + name
                    if len(new_name) <= IFORT_LIMIT:
                        break
                    name = self.__long_symbol_callback(name)
                else:
                    raise Exception("Too many attempts to rename long symbol")

                self.__renames[match.group('name').lower()] \
                    = self.__name.lower() + match.group('name').lower()

                new_line = match.group('head') + new_name \
                    + match.group('tail') + '\n'
                new_text.append(new_line)
                continue

            new_text.append(line)
        return ''.join(new_text)

    def rename_module_files(self, filename: Path) -> None:
        """
        Adds namespace to filename.

        This interacts with version control.

        If the file name already has the namespace it is skipped.

        :param filename: File to be renamed.
        """
        if filename.name.startswith(self.__name):
            return

        print("Renaming: " + str(filename))
        leafname = self.__name + filename.name
        new_filename = filename.with_name(leafname)
        run(['fcm', 'mv', str(filename), str(new_filename)], check=True)

        new_text = self.rename_module(new_filename.read_text())

        working_filename = new_filename.with_suffix('.wrk')
        working_filename.write_text(new_text)

        new_filename.unlink()
        working_filename.rename(new_filename)

    def write_renames(self, filename: Path):
        """
        Writes symbol mappings to file.

        :param filename: File to create.
        """
        with filename.open('wt') as file_handle:
            for key, value in self.__renames.items():
                print(f"{key} -> {value}", file=file_handle)


class Updater(ABC):
    """
    Parent for classes which update source files.
    """
    def __init__(self, map_file: Path):
        """
        :param map_file: File containing mappings.
        """
        self.__renames: Dict[str, str] = {}
        for line in map_file.read_text().splitlines():
            if line.strip() == '':  # Trap empty lines
                continue
            orig, new = re.split(r'\s*->\s*', line, maxsplit=1)
            self.__renames[orig] = new

    def rename(self, name: str) -> str:
        """
        Coverts symbol name to new name.
        """
        return self.__renames.get(name.lower()) or name

    def update_file(self, filename: Path) -> None:
        """
        Manages update of file.
        """
        if filename.is_symlink():
            return

        print(f"Updating ({self.__class__.__name__}): " + str(filename))
        new_text = self.update(filename.read_text())

        working_filename = filename.with_suffix('.wrk')
        working_filename.write_text(''.join(new_text))

        filename.unlink()
        working_filename.rename(filename)

    @abstractmethod
    def update(self, source: str) -> str:
        """
        Performs update of source file content.
        """
        raise NotImplementedError()


class UpdateUses(Updater):
    """
    Rewrites "use" lines with new module names.
    """
    def __init__(self, rename_file: Path):
        """
        Constructs a usage rewriter from rename mappings.

        :param rename_file: File holding symbol mappings.
        """
        super().__init__(rename_file)

        # Compile regular expressions for later.
        #
        self.__use_pattern = re_compile(
            r'(?P<indent>\s*)(?P<use>use\s+)(?P<module>\w+)'
            r'((?P<comma>\s*,)(?P<only_offset>\s*)(?P<only>only\s*:\s*)'
            r'(?P<only_list>(,?\s*[^&]+\s*)+))?(?P<continue>&.*)?',
            re.IGNORECASE
        )
        self.__continue_pattern = re_compile(
            r'(?P<indent>\s*)(?P<tail>.*)(?P<continue>&.*)?',
            re.IGNORECASE
        )

    def update(self, source: str) -> str:  #pylint: disable=too-many-branches
        indent = 0
        new_text: List[str] = []
        state = 'scanning'
        for line in source.splitlines(keepends=True):
            if state == 'scanning':
                match = self.__use_pattern.match(line)
                if match:
                    groups = match.groupdict()

                    new_name = self.rename(groups['module'])
                    new_line = groups['indent'] + groups['use']
                    new_line += new_name
                    if groups['only'] is not None:
                        offset_len = max(1,
                                         len(groups['only_offset'])
                                         + len(groups['module'])
                                         - len(new_name))
                        new_line += groups['comma'] \
                            + ' ' * offset_len \
                            + groups['only']
                        if groups['continue'] is not None:
                            indent = len(new_line)
                        if groups['only_list'] is not None:
                            new_line += groups['only_list']
                    if groups['continue'] is not None:
                        new_line += groups['continue']
                        state = 'continuing'

                    if not new_line.endswith('\n'):
                        new_line += '\n'
                    new_text.append(new_line)
                else:  # No match
                    new_text.append(line)
            elif state == 'continuing':
                match = self.__continue_pattern.match(line)
                if match:
                    groups = match.groupdict()

                    new_line = ' ' * indent
                    new_line += groups['tail']

                    if groups['continue'] is not None:
                        new_line += groups['continue']
                    else:
                        state = 'scanning'

                    if not new_line.endswith('\n'):
                        new_line += '\n'
                    new_text.append(new_line)
                else:  # No match
                    new_text.append(line)
                    state = 'scanning'
        return ''.join(new_text)
