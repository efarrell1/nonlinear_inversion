
      module run_star_extras

      use star_lib
      use star_def
      use const_def
      use crlibm_lib

      implicit none

      contains

        subroutine extras_controls(id, ierr)
           integer, intent(in) :: id
           integer, intent(out) :: ierr
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! this is the place to set any procedure pointers you want to change
           ! e.g., other_wind, other_mixing, other_energy  (see star_data.inc)

           ! Uncomment these lines if you wish to use the functions in this file,
           ! otherwise we use a null_ version which does nothing.
           ! s% extras_startup => extras_startup
           ! s% extras_start_step => extras_start_step
           ! s% extras_check_model => extras_check_model
           s% extras_finish_step => extras_finish_step
           ! s% extras_after_evolve => extras_after_evolve
           ! s% how_many_extra_history_columns => how_many_extra_history_columns
           ! s% data_for_extra_history_columns => data_for_extra_history_columns
           ! s% how_many_extra_profile_columns => how_many_extra_profile_columns
           ! s% data_for_extra_profile_columns => data_for_extra_profile_columns

           ! s% how_many_extra_history_header_items => how_many_extra_history_header_items
           ! s% data_for_extra_history_header_items => data_for_extra_history_header_items
           ! s% how_many_extra_profile_header_items => how_many_extra_profile_header_items
           ! s% data_for_extra_profile_header_items => data_for_extra_profile_header_items

           ! Once you have set the function pointers you want,
           ! then uncomment this (or set it in your star_job inlist)
           ! to disable the printed warning message,
            s% job% warn_run_star_extras =.false.

        end subroutine extras_controls


        integer function extras_finish_step(id, id_extra)
           integer, intent(in) :: id, id_extra
           integer :: ierr
           integer :: k
           real :: epsgrav
           real :: power_photo
           real :: diff_lum
           real :: max_eps_grav
           real :: max_eps_nuc
           real :: Lum
           real :: diff_lum_arr(10000)
           integer :: iphoto
           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return
           extras_finish_step = keep_going
           call store_extra_info(s)

           ! How to tell if a star is in thermal equilibrium
           if (s% x_logical_ctrl(1)) then

             iphoto = 10
             epsgrav = dot_product(s% dm(1:s% nz), s% eps_grav(1:s% nz))/3.8418d33

             power_photo = dot_product(s% dm(1:s% nz), s% eps_nuc_categories(iphoto,1:s% nz))/3.8418d33
             Lum = s% power_nuc_burn - power_photo

             max_eps_grav = maxval(abs(s% eps_grav))
             max_eps_nuc = maxval(abs(s% eps_nuc))

             do k = 1, s% nz - 2
               diff_lum_arr(k) = maxval((s% L(k + 1 : s% nz) / s% L(k)))
             end do


             write(*,*) abs(epsgrav / Lum), max_eps_grav/max_eps_nuc, maxval(diff_lum_arr)

             if (abs(epsgrav / Lum) < 0.01 .AND. max_eps_grav < max_eps_nuc * 0.1 .AND. maxval(diff_lum_arr) < 1.25) then
                 write(*,*) abs(epsgrav / Lum)
                 extras_finish_step = terminate
             end if
           end if

           if (extras_finish_step == terminate) s% termination_code = t_extras_finish_step

        end function extras_finish_step







        ! routines for saving and restoring extra data so can do restarts

           ! put these defs at the top and delete from the following routines
           !integer, parameter :: extra_info_alloc = 1
           !integer, parameter :: extra_info_get = 2
           !integer, parameter :: extra_info_put = 3


        subroutine alloc_extra_info(s)
           integer, parameter :: extra_info_alloc = 1
           type (star_info), pointer :: s
           call move_extra_info(s,extra_info_alloc)
        end subroutine alloc_extra_info


        subroutine unpack_extra_info(s)
           integer, parameter :: extra_info_get = 2
           type (star_info), pointer :: s
           call move_extra_info(s,extra_info_get)
        end subroutine unpack_extra_info


        subroutine store_extra_info(s)
           integer, parameter :: extra_info_put = 3
           type (star_info), pointer :: s
           call move_extra_info(s,extra_info_put)
        end subroutine store_extra_info


        subroutine move_extra_info(s,op)
           integer, parameter :: extra_info_alloc = 1
           integer, parameter :: extra_info_get = 2
           integer, parameter :: extra_info_put = 3
           type (star_info), pointer :: s
           integer, intent(in) :: op

           integer :: i, j, num_ints, num_dbls, ierr

           i = 0
           ! call move_int or move_flg
           num_ints = i

           i = 0
           ! call move_dbl

           num_dbls = i

           if (op /= extra_info_alloc) return
           if (num_ints == 0 .and. num_dbls == 0) return

           ierr = 0
           call star_alloc_extras(s% id, num_ints, num_dbls, ierr)
           if (ierr /= 0) then
              write(*,*) 'failed in star_alloc_extras'
              write(*,*) 'alloc_extras num_ints', num_ints
              write(*,*) 'alloc_extras num_dbls', num_dbls
              stop 1
           end if

           contains

           subroutine move_dbl(dbl)
              real(dp) :: dbl
              i = i+1
              select case (op)
              case (extra_info_get)
                 dbl = s% extra_work(i)
              case (extra_info_put)
                 s% extra_work(i) = dbl
              end select
           end subroutine move_dbl

           subroutine move_int(int)
              integer :: int
              i = i+1
              select case (op)
              case (extra_info_get)
                 int = s% extra_iwork(i)
              case (extra_info_put)
                 s% extra_iwork(i) = int
              end select
           end subroutine move_int

           subroutine move_flg(flg)
              logical :: flg
              i = i+1
              select case (op)
              case (extra_info_get)
                 flg = (s% extra_iwork(i) /= 0)
              case (extra_info_put)
                 if (flg) then
                    s% extra_iwork(i) = 1
                 else
                    s% extra_iwork(i) = 0
                 end if
              end select
           end subroutine move_flg

        end subroutine move_extra_info


      end module run_star_extras
