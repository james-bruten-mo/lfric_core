##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import subprocess
from pathlib import Path
from subprocess import run
from textwrap import dedent

from pytest import fixture, raises

from ..modules.namespacerator import AddNamespace, UpdateUses


class TestRenamer:
    def test_rename_module(self, tmp_path: Path):
        source = dedent(
            """
            ! First module is vanilla
            !
            module first_mod
            end module first_mod
            !
            ! Second module has old school style.
            !
            MODULE second_mod
            ENDMODULE second_mod
            !
            ! Third module has missing name on "end" statement.
            !
            module third_mod
            end module
            """
        )
        expected = dedent(
            """
            ! First module is vanilla
            !
            module foo_first_mod
            end module foo_first_mod
            !
            ! Second module has old school style.
            !
            MODULE foo_second_mod
            ENDMODULE foo_second_mod
            !
            ! Third module has missing name on "end" statement.
            !
            module foo_third_mod
            end module
            """
        )

        test_unit = AddNamespace("foo")
        result = test_unit.rename_module(source)
        assert result == expected

        repo_path = tmp_path / 'repo'
        run(['svnadmin', 'create', repo_path])
        working_copy_path = tmp_path / 'wc'
        run(['svn', 'checkout', f'file://{str(repo_path)}', working_copy_path])

        source_file = working_copy_path / 'source.f90'
        source_file.write_text(source)
        run(['svn', 'add', source_file], cwd=working_copy_path)
        run(['svn', 'commit', '-m', 'Test file added'], cwd=working_copy_path)

        test_unit.rename_module_files(source_file)
        process = run(['svn', 'status'], cwd=working_copy_path,
                      stdout=subprocess.PIPE, encoding='utf-8')
        assert process.stdout == dedent(
            '''
            A  +    foo_source.f90
                    > moved from source.f90
            D       source.f90
                    > moved to foo_source.f90
            '''
        ).lstrip()

        expected_filename = working_copy_path / 'foo_source.f90'
        assert expected_filename.read_text(encoding='utf8') == expected

        map_path = tmp_path / 'map.txt'
        test_unit.write_renames(map_path)
        assert map_path.read_text(encoding='utf8') == dedent(
            '''
            first_mod -> foo_first_mod
            second_mod -> foo_second_mod
            third_mod -> foo_third_mod
            '''
        ).lstrip()

    def test_module_procedure(self):
        source = dedent(
            """
            module with_interface_mod
              interface foo_if
                module procedure bar
                module procedure baz
              end interface foo_if
            end module with_interface_mod
            """
        )
        test_unit = AddNamespace('cheese')
        result = test_unit.rename_module(source)
        assert result == dedent(
            """
            module cheese_with_interface_mod
              interface foo_if
                module procedure bar
                module procedure baz
              end interface foo_if
            end module cheese_with_interface_mod
            """
        )

    @fixture(scope='class')
    def long_name_module_source(self) -> str:
        return dedent(
            """
        module very_very_extremely_ludicrously_ridiculously_long_name_mod
        end module very_very_extremely_ludicrously_ridiculously_long_name_mod
        """
        )

    def test_default_ifort_length(self, long_name_module_source):
        # Default behaviour is to raise an exception.
        #
        test_unit = AddNamespace('bfort')
        with raises(Exception):
            _ = test_unit.rename_module(long_name_module_source)

    def test_override_ifort_length(self, long_name_module_source):
        # Override default behaviour with callback.
        #
        test_unit = AddNamespace('ifort', long_symbol_callback=lambda x: x[5:])
        result = test_unit.rename_module(long_name_module_source)
        assert result == dedent(
            """
        module ifort_very_extremely_ludicrously_ridiculously_long_name_mod
        end module ifort_very_extremely_ludicrously_ridiculously_long_name_mod
        """
        )

    def test_give_up_ifort_length(self, long_name_module_source):
        # Override default behaviour with faulty callback.
        #
        test_unit = AddNamespace('cfort', long_symbol_callback=lambda x: x)
        with raises(Exception):
            _ = test_unit.rename_module(long_name_module_source)


class TestUpdateUsage:
    def test_update_usage(self, tmp_path: Path):
        rename_file = tmp_path / 'map.txt'
        rename_file.write_text(dedent(
            """
            simple_mod  -> bar_simple_mod
            specific_mod->bar_specific_mod
            big_mod     -> bar_big_mod
            bigerer_mod -> bar_bigerer_mod
            understandable -> bar_understandable
            thing_the_second -> bar_thing_the_second
            medium_mod  -> bar_medium_mod
            """
        ))
        source = dedent(
            """
            ! Ersatz copyright statement.
            !
            module test_mod
                ! Use some things which weren't renamed
                use not_this_mod
                use nor_this_mod, only : unchanged_thing
                use never_this, only : thing_one, thing_two
                ! Use some things which were renamed
                use simple_mod
                use specific_mod, only : changed_thing
                use big_mod,      only: thing_the_first,   &
                                        thing_the_second,  &
                                        thing_the_third
                use bigerer_mod, only: well_named => argle_bargle, &
                                       understandable => biggly
            contains
                subroutine stuff
                    integer :: used_colours(4)
                    used_colours = (/1, 2, 3, 4/)
                end subroutine
            end module test_mod
            """
        )
        expected = dedent(
            """
            ! Ersatz copyright statement.
            !
            module test_mod
                ! Use some things which weren't renamed
                use not_this_mod
                use nor_this_mod, only : unchanged_thing
                use never_this, only : thing_one, thing_two
                ! Use some things which were renamed
                use bar_simple_mod
                use bar_specific_mod, only : changed_thing
                use bar_big_mod,  only: thing_the_first,   &
                                        thing_the_second,  &
                                        thing_the_third
                use bar_bigerer_mod, only: well_named => argle_bargle, &
                                           understandable => biggly
            contains
                subroutine stuff
                    integer :: used_colours(4)
                    used_colours = (/1, 2, 3, 4/)
                end subroutine
            end module test_mod
            """
        )
        test_unit = UpdateUses(rename_file)
        result = test_unit.update(source)
        assert result == expected

        source_file = tmp_path / 'source.f90'
        source_file.write_text(source)
        test_unit.update_file(source_file)
        assert source_file.read_text() == expected
