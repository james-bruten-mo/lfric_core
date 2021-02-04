!!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!
!> @brief ENUM containing different field synonyms methods
!> @details This enumerator contains all the possible field synonyms
!> methods for an LFRic field. It is used by meta_data_mod.f90 for
!> specifying the recommended field synonyms type

module field_synonyms_enum_mod

  implicit none
  !> THE NUMBERS IN THIS FILE ARE ARBITARY.
  !> DO NOT ATTEMPT TO USE THEM IN ANY CODE
  enum, bind(c)

    enumerator :: AMIP = 1001
    enumerator :: CF = 1002
    enumerator :: CMIP6 = 1003
    enumerator :: GRIB = 1004
    enumerator :: STASH = 1005

  end enum

end module field_synonyms_enum_mod