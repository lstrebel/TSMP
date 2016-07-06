! This code is not necessarily under the DART copyright ...
!
! DART $Id: $

module model_mod

! This module provides routines to work with COSMO data
! files in the DART framework
!
! Based on a version by : Jan D. Keller 2011-09-15
!         Meteorological Institute, University of Bonn, Germany

use        types_mod, only : i4, r4, r8, digits12, SECPERDAY, MISSING_R8,          &
                             rad2deg, deg2rad, PI, obstypelength

use time_manager_mod, only : time_type, set_time, set_date, get_date, get_time,&
                             print_time, print_date, set_calendar_type,        &
                             operator(*),  operator(+), operator(-),           &
                             operator(>),  operator(<), operator(/),           &
                             operator(/=), operator(<=)

use   cosmo_data_mod, only : cosmo_meta, cosmo_hcoord, cosmo_non_state_data,   &
                             get_cosmo_info, get_data_from_binary,             &
                             set_vertical_coords, grib_header_type,            &
                             model_dims, record_length

use     location_mod, only : location_type, get_dist, query_location,          &
                             get_close_maxdist_init, get_close_type,           &
                             set_location, get_location, horiz_dist_only,      &
                             vert_is_undef,        VERTISUNDEF,                &
                             vert_is_surface,      VERTISSURFACE,              &
                             vert_is_level,        VERTISLEVEL,                &
                             vert_is_pressure,     VERTISPRESSURE,             &
                             vert_is_height,       VERTISHEIGHT,               &
                             vert_is_scale_height, VERTISSCALEHEIGHT,          &
                             get_close_obs_init, get_close_obs

use    utilities_mod, only : register_module, error_handler,                   &
                             E_ERR, E_WARN, E_MSG, logfileunit, get_unit,      &
                             nc_check, do_output, to_upper,                    &
                             find_namelist_in_file, check_namelist_read,       &
                             open_file, close_file, file_exist,                &
                             find_textfile_dims, file_to_text,                 &
                             do_nml_file, do_nml_term

use     obs_kind_mod, only : KIND_U_WIND_COMPONENT,       &
                             KIND_V_WIND_COMPONENT,       &
                             KIND_VERTICAL_VELOCITY,      &
                             KIND_TEMPERATURE,            &
                             KIND_PRESSURE,               &
                             KIND_PRESSURE_PERTURBATION,  &
                             KIND_SPECIFIC_HUMIDITY,      &
                             KIND_CLOUD_LIQUID_WATER,     &
                             KIND_CLOUD_ICE,              &
                             KIND_SURFACE_ELEVATION,      &
                             KIND_SURFACE_GEOPOTENTIAL,   &
                             paramname_length,            &
                             get_raw_obs_kind_index,      &
                             get_raw_obs_kind_name

use    random_seq_mod, only: random_seq_type, init_random_seq, random_gaussian

use byte_mod, only: to_float1,from_float1,word_to_byte,byte_to_word_signed,concat_bytes1

use netcdf

implicit none
private

! version controlled file description for error handling, do not edit
character(len=256), parameter :: source   = "$URL: model_mod.f90 $"
character(len=32 ), parameter :: revision = "$Revision: none $"
character(len=128), parameter :: revdate  = "$Date: none $"

character(len=256) :: string1, string2, string3
logical, save :: module_initialized = .false.

! these routines must be public and you cannot change
! the arguments - they will be called *from* the DART code.

public :: get_model_size,         &
          adv_1step,              &
          get_state_meta_data,    &
          model_interpolate,      &
          get_model_time_step,    &
          static_init_model,      &
          end_model,              &
          init_time,              &
          init_conditions,        &
          nc_write_model_atts,    &
          nc_write_model_vars,    &
          pert_model_state,       &
          get_close_maxdist_init, &
          get_close_obs_init,     &
          get_close_obs,          &
          ens_mean_for_model

! generally useful routines for various support purposes.
! the interfaces here can be changed as appropriate.

!  public  :: grib_to_sv
!  public  :: sv_to_grib

public :: get_state_time,     &
          get_state_vector,   &
          write_grib_file,    &
          get_cosmo_filename, &
          write_state_times

INTERFACE sv_to_field
   MODULE PROCEDURE sv_to_field_2d
   MODULE PROCEDURE sv_to_field_3d
END INTERFACE

! TODO FIXME  ultimately the dart_variable_info type will be removed
type dart_variable_info
   character(len=16)    :: varname_short
   character(len=256)   :: varname_long
   character(len=32)    :: units
   logical              :: is_present
   integer              :: nx
   integer              :: ny
   integer              :: nz
   real(r8),allocatable :: vertical_level(:)
   integer              :: vertical_coordinate
   integer              :: horizontal_coordinate
   integer,allocatable  :: state_vector_sindex(:) ! starting index in state vector for every vertical level
   integer,allocatable  :: cosmo_state_index(:)   ! index in cosmo state of every vertical level
end type dart_variable_info

! Codes for restricting the range of a variable
integer, parameter :: BOUNDED_NONE  = 0 ! ... unlimited range
integer, parameter :: BOUNDED_BELOW = 1 ! ... minimum, but no maximum
integer, parameter :: BOUNDED_ABOVE = 2 ! ... maximum, but no minimum
integer, parameter :: BOUNDED_BOTH  = 3 ! ... minimum and maximum

integer :: nfields
integer, parameter :: max_state_variables = 80
integer, parameter :: num_state_table_columns = 7
character(len=obstypelength) :: variable_table(max_state_variables, num_state_table_columns)

! Codes for interpreting the columns of the variable_table
integer, parameter :: VT_GRIBVERSIONINDX = 1 ! ... the version of the grib table being used
integer, parameter :: VT_GRIBVARINDX     = 2 ! ... variable name
integer, parameter :: VT_VARNAMEINDX     = 3 ! ... netcdf variable name
integer, parameter :: VT_KINDINDX        = 4 ! ... DART kind
integer, parameter :: VT_MINVALINDX      = 5 ! ... minimum value if any
integer, parameter :: VT_MAXVALINDX      = 6 ! ... maximum value if any
integer, parameter :: VT_STATEINDX       = 7 ! ... update (state) or not

! Everything needed to describe a variable
!>@ TODO FIXME remove the unused netcdf bits ... we're working with binary only

type progvartype
   private
   character(len=NF90_MAX_NAME) :: varname
   character(len=NF90_MAX_NAME) :: long_name
   character(len=NF90_MAX_NAME) :: units
   character(len=obstypelength) :: dimnames(NF90_MAX_VAR_DIMS)
   integer  :: tableID     ! grib table version
   integer  :: variableID  ! variable ID in grib table
   integer  :: levtypeID   ! the kind of vertical coordinate system
   integer  :: dimlens(NF90_MAX_VAR_DIMS)
   integer  :: numdims
   integer  :: maxlevels
   integer  :: varsize     ! prod(dimlens(1:numdims))
   integer  :: index1      ! location in dart state vector of first occurrence
   integer  :: indexN      ! location in dart state vector of last  occurrence
   integer  :: dart_kind
   integer  :: rangeRestricted
   real(r8) :: minvalue
   real(r8) :: maxvalue
   character(len=paramname_length) :: kind_string
   logical  :: update
end type progvartype

type(progvartype), dimension(max_state_variables) :: progvar

! Dimensions (from the netCDF file) specifying the grid shape, etc.
integer :: ntime
integer :: nbnds
integer :: nrlon
integer :: nrlat
integer :: nsrlon
integer :: nsrlat
integer :: nlevel
integer :: nlevel1
integer :: nsoil
integer :: nsoil1

!> @TODO FIXME ... check to make sure the storage order in netCDF is the
!> same as the storage order of the binary variables.
!> vcoord has units, etc. and will require a transformation

real(r8), allocatable ::   lon(:,:) ! longitude degrees_east
real(r8), allocatable ::   lat(:,:) ! latitude  degrees_north
real(r8), allocatable :: slonu(:,:) ! staggered U-wind longitude degrees_east
real(r8), allocatable :: slatu(:,:) ! staggered U-wind latitude  degrees_north
real(r8), allocatable :: slonv(:,:) ! staggered V-wind longitude degrees_east
real(r8), allocatable :: slatv(:,:) ! staggered V-wind latitude  degrees_north

! just mimic what is in the netCDF variable
!       float vcoord(level1) ;
!               vcoord:long_name = "Height-based hybrid Gal-Chen coordinate" ;
!               vcoord:units = "Pa" ;
!               vcoord:ivctype = 2 ;
!               vcoord:irefatm = 2 ;
!               vcoord:p0sl = 100000. ;
!               vcoord:t0sl = 300. ;
!               vcoord:dt0lp = 42. ;
!               vcoord:vcflat = 11000. ;
!               vcoord:delta_t = 75. ;
!               vcoord:h_scal = 10000. ;

type verticalobject
   private
   real(r8), allocatable :: level1(:)
   character(len=128)    :: long_name = "Height-based hybrid Gal-Chen coordinate"
   character(len=32)     :: units = "Pa"
   integer               :: ivctype = 2
   integer               :: irefatm = 2
   real(r8)              :: p0sl = 100000.
   real(r8)              :: t0sl = 300.
   real(r8)              :: dt0lp = 42.
   real(r8)              :: vcflat = 11000.
   real(r8)              :: delta_t = 75.
   real(r8)              :: h_scal = 10000.
end type verticalobject

type(verticalobject) :: vcoord

! track which variables to update
! As the binary file gets read, we need to compare the 'slab' to see if
! it is a slab we want to update.
!>@ TODO FIXME after we populate progvar, ordered_update_list needs
!> to be populated with the order to stride through progvar to
!> write to the output file.

type update_items
   private
   integer :: tableID
   integer :: variableID
   integer :: levtype
   integer :: nlevels
end type update_items
type(update_items) :: ordered_update_list(max_state_variables)

integer, parameter             :: n_state_vector_vars=8
integer, parameter             :: n_non_state_vars=1

type(cosmo_meta),allocatable   :: cosmo_slabs(:)
type(cosmo_hcoord)             :: cosmo_lonlat(3) ! 3 is for the stagger
integer                        :: nslabs

! things which can/should be in the model_nml

character(len=256) :: cosmo_restart_file           = "cosmo_restart_file"
character(len=256) :: cosmo_netcdf_file            = "cosmo_netcdf_file"
integer            :: assimilation_period_days     = 0
integer            :: assimilation_period_seconds  = 60
integer            :: model_dt                     = 40
logical            :: output_1D_state_vector       = .FALSE.
real(r8)           :: model_perturbation_amplitude = 0.1
integer            :: debug                        = 0
character(len=obstypelength) :: variables(max_state_variables*num_state_table_columns) = ' '

namelist /model_nml/             &
   cosmo_restart_file,           &
   cosmo_netcdf_file,            &
   assimilation_period_days,     &
   assimilation_period_seconds,  &
   model_dt,                     &
   model_perturbation_amplitude, &
   output_1D_state_vector,       &
   model_dims,                   &
   record_length,                &
   debug,                        &
   variables

integer                        :: model_size
type(time_type)                :: model_timestep ! smallest time to adv model

integer, parameter             :: n_max_kinds=400

integer                        :: allowed_state_vector_vars(n_state_vector_vars)
integer                        :: allowed_non_state_vars(1:n_max_kinds)
logical                        :: is_allowed_state_vector_var(n_max_kinds)
logical                        :: is_allowed_non_state_var(n_max_kinds)
type(dart_variable_info)       :: state_vector_vars(1:n_max_kinds)
type(cosmo_non_state_data)     :: non_state_data

real(r8),allocatable           :: state_vector(:)

real(r8),allocatable           :: ens_mean(:)

type(random_seq_type) :: random_seq

type(time_type)                :: cosmo_fc_time
type(time_type)                :: cosmo_an_time

type(grib_header_type),allocatable  :: grib_header(:)

contains


!------------------------------------------------------------------------
!>


function get_model_size()

  integer :: get_model_size

  call error_handler(E_ERR,'get_model_size','routine not written',source,revision,revdate)

  if ( .not. module_initialized ) call static_init_model

  get_model_size = model_size

end function get_model_size


!------------------------------------------------------------------------
!> Called to do one-time initialization of the model.
!>
!> All the grid information comes from the COSMOS netCDF file
!>@ TODO FIXME All the variable information comes from all over the place
!> Not actually reading in the state, that is done in get_state_vector()

subroutine static_init_model()

! Local variables - all the important ones have module scope

integer               :: ivar, index1, indexN
integer               :: iunit, io, islab, ikind, sv_length
real(r8), allocatable :: datmat(:,:)
real(r8), parameter   :: g = 9.80665_r8

if ( module_initialized ) return ! only need to do this once.

module_initialized = .TRUE.

! read the DART namelist for this model
call find_namelist_in_file('input.nml', 'model_nml', iunit)
read(iunit, nml = model_nml, iostat = io)
call check_namelist_read(iunit, io, 'model_nml')

! Record the namelist values used for the run
if (do_nml_file()) write(logfileunit, nml=model_nml)
if (do_nml_term()) write(     *     , nml=model_nml)

call set_calendar_type('Gregorian')

! Get the dimensions of the grid and the grid variables from the netCDF file.
call get_cosmo_grid(cosmo_netcdf_file)

! rectify the user input for the variables to include in the DART state vector
call parse_variable_table(variables, nfields, variable_table)

! do ivar = 1,nfields
do ivar = 1,1
!   call set_variable_attributes(ivar)
    call set_variable_binary_properties(ivar)
    if (debug > 5 .and. do_output()) call progvar_summary()
enddo

call error_handler(E_ERR,'static_init_model','routine not finished',source,revision,revdate)

call get_cosmo_info(cosmo_restart_file, cosmo_slabs, cosmo_lonlat, grib_header, &
                      is_allowed_state_vector_var, cosmo_fc_time)


call set_allowed_state_vector_vars()

  state_vector_vars(:)%is_present=.false.
  non_state_data%orography_present=.false.
  non_state_data%pressure_perturbation_present=.false.

  model_size = maxval(cosmo_slabs(:)%dart_eindex)
  nslabs     = size(cosmo_slabs,1)

  sv_length = 0
  do islab = 1,nslabs
    ikind = cosmo_slabs(islab)%dart_kind
    if ( ikind > 0) then
      if (is_allowed_state_vector_var(ikind).OR.(ikind==KIND_PRESSURE_PERTURBATION)) then

        sv_length = sv_length + cosmo_slabs(islab)%dims(1)*cosmo_slabs(islab)%dims(2)

      endif
    endif
  enddo

  allocate(state_vector(1:sv_length))

  ! cycle through all GRIB records
  ! one record corresponds to one horizontal field
  do islab=1,nslabs
    ikind=cosmo_slabs(islab)%dart_kind

    ! check if variable is a possible state vector variable
    if (ikind>0) then
      if (is_allowed_state_vector_var(ikind)) then

        ! check if state vector variable information has not already been read
        ! e.g. for another vertical level
        if (.not. state_vector_vars(ikind)%is_present) then
          ! assign the variable information
          state_vector_vars(ikind)%is_present   = .true.
          state_vector_vars(ikind)%varname_short= cosmo_slabs(islab)%varname_short
          state_vector_vars(ikind)%varname_long = cosmo_slabs(islab)%varname_long
          state_vector_vars(ikind)%units        = cosmo_slabs(islab)%units
          state_vector_vars(ikind)%nx           = cosmo_slabs(islab)%dims(1)
          state_vector_vars(ikind)%ny           = cosmo_slabs(islab)%dims(2)
          state_vector_vars(ikind)%nz           = cosmo_slabs(islab)%dims(3)
          state_vector_vars(ikind)%horizontal_coordinate=cosmo_slabs(islab)%hcoord_type
          if (state_vector_vars(ikind)%nz>1) then
            state_vector_vars(ikind)%vertical_coordinate=VERTISLEVEL
          else
            state_vector_vars(ikind)%vertical_coordinate=VERTISSURFACE
          endif

          allocate(state_vector_vars(ikind)%vertical_level(     1:state_vector_vars(ikind)%nz))
          allocate(state_vector_vars(ikind)%state_vector_sindex(1:state_vector_vars(ikind)%nz))
          allocate(state_vector_vars(ikind)%cosmo_state_index(  1:state_vector_vars(ikind)%nz))

        endif

        ! set vertical information for this vertical level (record/slab)
        state_vector_vars(ikind)%vertical_level(     cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_level
        state_vector_vars(ikind)%state_vector_sindex(cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_sindex
        state_vector_vars(ikind)%cosmo_state_index(  cosmo_slabs(islab)%ilevel)=islab

      endif

    ! check for non state vector datmat (e.g. surface elevation) needed to run DART
      if (is_allowed_non_state_var(ikind)) then
        if (ikind==KIND_SURFACE_ELEVATION) then
          allocate(datmat(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
          datmat=get_data_from_binary(cosmo_restart_file,grib_header(islab),cosmo_slabs(islab)%dims(1),cosmo_slabs(islab)%dims(2))
          if (.not. allocated(non_state_data%surface_orography)) then
            allocate(non_state_data%surface_orography(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
          endif
          non_state_data%surface_orography(:,:)=datmat(:,:)
          deallocate(datmat)
          non_state_data%orography_present=.true.
        endif
        if ((ikind==KIND_SURFACE_GEOPOTENTIAL).and.(.not. allocated(non_state_data%surface_orography))) then
          allocate(datmat(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
          datmat=get_data_from_binary(cosmo_restart_file,grib_header(islab),cosmo_slabs(islab)%dims(1),cosmo_slabs(islab)%dims(2))
          allocate(non_state_data%surface_orography(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
          non_state_data%surface_orography(:,:)=datmat(:,:)/g
          deallocate(datmat)
          non_state_data%orography_present=.true.
        endif
        if (ikind==KIND_PRESSURE_PERTURBATION) then
          allocate(datmat(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2)))
          datmat=get_data_from_binary(cosmo_restart_file,grib_header(islab),cosmo_slabs(islab)%dims(1),cosmo_slabs(islab)%dims(2))
          if (.not. allocated(non_state_data%pressure_perturbation)) then
            allocate(non_state_data%pressure_perturbation(1:cosmo_slabs(islab)%dims(1),1:cosmo_slabs(islab)%dims(2),1:cosmo_slabs(islab)%dims(3)))
          endif
          non_state_data%pressure_perturbation(:,:,cosmo_slabs(islab)%ilevel)=datmat(:,:)
          deallocate(datmat)

          if (.not. state_vector_vars(KIND_PRESSURE)%is_present) then
            ! assign the pressure variable information
            state_vector_vars(KIND_PRESSURE)%is_present   = .true.
            state_vector_vars(KIND_PRESSURE)%varname_short= cosmo_slabs(islab)%varname_short
            state_vector_vars(KIND_PRESSURE)%varname_long = cosmo_slabs(islab)%varname_long
            state_vector_vars(KIND_PRESSURE)%units        = cosmo_slabs(islab)%units
            state_vector_vars(KIND_PRESSURE)%nx           = cosmo_slabs(islab)%dims(1)
            state_vector_vars(KIND_PRESSURE)%ny           = cosmo_slabs(islab)%dims(2)
            state_vector_vars(KIND_PRESSURE)%nz           = cosmo_slabs(islab)%dims(3)
            state_vector_vars(KIND_PRESSURE)%horizontal_coordinate=cosmo_slabs(islab)%hcoord_type
            state_vector_vars(KIND_PRESSURE)%vertical_coordinate=VERTISLEVEL
            allocate(state_vector_vars(KIND_PRESSURE)%vertical_level(     1:state_vector_vars(KIND_PRESSURE)%nz))
            allocate(state_vector_vars(KIND_PRESSURE)%state_vector_sindex(1:state_vector_vars(KIND_PRESSURE)%nz))
            allocate(state_vector_vars(KIND_PRESSURE)%cosmo_state_index(  1:state_vector_vars(KIND_PRESSURE)%nz))
          endif

          ! set vertical information for this vertical level (record/slab)
          state_vector_vars(KIND_PRESSURE)%vertical_level(     cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_level
          state_vector_vars(KIND_PRESSURE)%state_vector_sindex(cosmo_slabs(islab)%ilevel)=cosmo_slabs(islab)%dart_sindex
          state_vector_vars(KIND_PRESSURE)%cosmo_state_index(  cosmo_slabs(islab)%ilevel)=islab

          non_state_data%pressure_perturbation_present=.true.
        endif

      endif
    endif
  enddo

  ! set up the vertical coordinate system information
  !   search for one 3D variable (U-wind component should be contained in every analysis file)
  setlevel : do islab=1,nslabs
    if (cosmo_slabs(islab)%dart_kind==KIND_U_WIND_COMPONENT) then

      ! calculate the vertical coordinates for every grid point
      call set_vertical_coords(grib_header(islab),non_state_data,state_vector_vars(KIND_PRESSURE)%state_vector_sindex(:),state_vector)

      exit setlevel
    endif
  enddo setlevel

  return
end subroutine static_init_model


!------------------------------------------------------------------------
!>


subroutine get_state_meta_data(index_in, location, var_type)

integer, intent(in)            :: index_in
type(location_type)            :: location
integer, optional, intent(out) :: var_type

integer   :: islab,var,hindex,dims(3)
real(r8)  :: mylon,mylat,vloc

call error_handler(E_ERR,'get_state_meta_data','routine not written',source,revision,revdate)

if (.NOT. module_initialized) CALL static_init_model()

var = -1

findindex : DO islab=1,nslabs
  IF ((index_in >= cosmo_slabs(islab)%dart_sindex) .AND. (index_in <= cosmo_slabs(islab)%dart_eindex)) THEN
    var      = islab
    hindex   = index_in - cosmo_slabs(islab)%dart_sindex + 1
    var_type = cosmo_slabs(islab)%dart_kind
    dims     = cosmo_slabs(islab)%dims
    vloc     = cosmo_slabs(islab)%dart_level
    mylon    = cosmo_lonlat(cosmo_slabs(islab)%hcoord_type)%lon(hindex)
    mylat    = cosmo_lonlat(cosmo_slabs(islab)%hcoord_type)%lat(hindex)
    location = set_location(mylon,mylat,vloc,VERTISLEVEL)
    EXIT findindex
  endif
enddo findindex

IF( var == -1 ) THEN
  write(string1,*) 'Problem, cannot find base_offset, index_in is: ', index_in
  call error_handler(E_ERR,'get_state_meta_data',string1,source,revision,revdate)
ENDIF

end subroutine get_state_meta_data


!------------------------------------------------------------------------
!>


  function get_model_time_step()
  ! Returns the smallest increment of time that we want to advance the model.
  ! This defines the minimum assimilation interval.
  ! It is NOT the dynamical timestep of the model.

    type(time_type) :: get_model_time_step

    call error_handler(E_ERR,'get_model_time_step','routine not written',source,revision,revdate)

    if ( .not. module_initialized ) call static_init_model

    model_timestep      = set_time(model_dt)
    get_model_time_step = model_timestep
    return

  end function get_model_time_step


!------------------------------------------------------------------------
!>


  subroutine model_interpolate(x, location, obs_type, interp_val, istatus)

  ! Error codes:
  ! istatus = 99 : unknown error
  ! istatus = 10 : observation type is not in state vector
  ! istatus = 15 : observation lies outside the model domain (horizontal)
  ! istatus = 16 : observation lies outside the model domain (vertical)
  ! istatus = 19 : observation vertical coordinate is not supported

  ! Passed variables

  real(r8),            intent(in)  :: x(:)
  type(location_type), intent(in)  :: location
  integer,             intent(in)  :: obs_type
  real(r8),            intent(out) :: interp_val
  integer,             intent(out) :: istatus

  ! Local storage

  real(r8)             :: point_coords(1:3)

  integer              :: i,j,hbox(2,2),n,vbound(2),sindex
  real(r8)             :: hbox_weight(2,2),hbox_val(2,2),hbox_lon(2,2),hbox_lat(2,2)
  real(r8)             :: vbound_weight(2),val1,val2

  call error_handler(E_ERR,'model_interpolate','routine not written',source,revision,revdate)

  IF ( .not. module_initialized ) call static_init_model

  interp_val = MISSING_R8     ! the DART bad value flag
  istatus = 99                ! unknown error

  ! FIXME ... want some sort of error message here?
  if ( .not. state_vector_vars(obs_type)%is_present) then
     istatus=10
     return
  endif

  ! horizontal interpolation

  n = size(cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon,1)

  point_coords(1:3) = get_location(location)

  ! Find grid indices of box enclosing the observation location
  call get_enclosing_grid_box_lonlat(cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon,&
                                     cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat,&
                                     point_coords(1:2),n,state_vector_vars(obs_type)%nx,                 &
                                     state_vector_vars(obs_type)%ny, hbox, hbox_weight)

  if (hbox(1,1)==-1) then
     istatus=15
     return
  endif

  ! determine vertical level above and below obsevation
  call get_vertical_boundaries(hbox, hbox_weight, obs_type, query_location(location,'which_vert'),&
                               point_coords(3), vbound, vbound_weight, istatus)

  ! check if observation is in vertical domain and vertical coordinate system is supported
  ! FIXME istatus value?
  if (vbound(1)==-1) then
     return
  endif

  ! Perform a bilinear interpolation from the grid box to the desired location
  ! for the level above and below the observation

  sindex=state_vector_vars(obs_type)%state_vector_sindex(vbound(1))

  do i=1,2
  do j=1,2
     hbox_val(i,j)=x(sindex+hbox(i,j)-1)
     hbox_lon(i,j)=cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lon(hbox(i,j))
     hbox_lat(i,j)=cosmo_lonlat(state_vector_vars(obs_type)%horizontal_coordinate)%lat(hbox(i,j))
  enddo
  enddo

  call bilinear_interpolation(hbox_val,hbox_lon,hbox_lat,point_coords,val1)

  sindex=state_vector_vars(obs_type)%state_vector_sindex(vbound(2))
  do i=1,2
  do j=1,2
     hbox_val(i,j)=x(sindex+hbox(i,j)-1)
  enddo
  enddo

  call bilinear_interpolation(hbox_val,hbox_lon,hbox_lat,point_coords,val2)

  ! vertical interpolation of horizontally interpolated values

  interp_val=val1*vbound_weight(1)+val2*vbound_weight(2)
  istatus=0

  end subroutine model_interpolate


!------------------------------------------------------------------------
!>
!> Returns a model state vector, x, that is some sort of appropriate
!> initial condition for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.

subroutine init_conditions(x)

real(r8), intent(out) :: x(:)

if ( .not. module_initialized ) call static_init_model

write(string1,*)'Cannot initialize COSMO time via subroutine call.'
write(string2,*)'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR, 'init_conditions', string1, &
           source, revision, revdate, text2=string2)

x = 0.0_r8  ! suppress compiler warnings about unused variables

end subroutine init_conditions


!------------------------------------------------------------------------
!> Companion interface to init_conditions. Returns a time that is somehow
!> appropriate for starting up a long integration of the model.
!> At present, this is only used if the namelist parameter
!> start_from_restart is set to .false. in the program perfect_model_obs.


subroutine init_time(time)

type(time_type), intent(out) :: time

if ( .not. module_initialized ) call static_init_model

write(string1,*)'Cannot initialize COSMO time via subroutine call.'
write(string2,*)'input.nml:start_from_restart cannot be FALSE'
call error_handler(E_ERR, 'init_time', string1, &
           source, revision, revdate, text2=string2)

time = set_time(0,0) ! suppress compiler warnings about unused variables

end subroutine init_time


!------------------------------------------------------------------------
!> As COSMO can only be advanced as a separate executable,
!> this is a NULL INTERFACE and is a fatal error if invoked.

subroutine adv_1step(x, time)

real(r8),        intent(inout) :: x(:)
type(time_type), intent(in)    :: time

if ( .not. module_initialized ) call static_init_model

if (do_output()) then
    call print_time(time,'NULL interface adv_1step (no advance) DART time is')
    call print_time(time,'NULL interface adv_1step (no advance) DART time is',logfileunit)
endif

write(string1,*) 'Cannot advance COSMO with a subroutine call; async cannot equal 0'
call error_handler(E_ERR,'adv_1step',string1,source,revision,revdate)

end subroutine adv_1step


!------------------------------------------------------------------------
!>


subroutine end_model()

deallocate(cosmo_slabs)
deallocate(state_vector)
deallocate(lon,lat,slonu,slatu,slonv,slatv)
deallocate(vcoord%level1)

end subroutine end_model


!------------------------------------------------------------------------
!>


  function nc_write_model_atts( ncFileID ) result (ierr)

    integer, intent(in)  :: ncFileID      ! netCDF file identifier
    integer              :: ierr          ! return value of function

    integer              :: nDimensions, nVariables, nAttributes, unlimitedDimID, TimeDimID
    integer              :: StateVarDimID   ! netCDF pointer to state variable dimension (model size)
    integer              :: MemberDimID     ! netCDF pointer to dimension of ensemble    (ens_size)
    integer              :: LineLenDimID
    integer              :: StateVarVarID,StateVarID,VarID
    integer              :: ikind,ndims,idim,dims(100),nx,ny,nz,i
    character(len=6)     :: ckind

    integer              :: lonDimID, latDimID, levDimID, wlevDimID
    integer              :: lonVarID, latVarID, ulonVarID, ulatVarID, vlonVarID, vlatVarID
    integer              :: levVarID, wlevVarID

    character(len=128)   :: filename
    real(r8)             :: levs(1:500),wlevs(1:501)
    real(r8),allocatable :: data2d(:,:)

    character(len=8)      :: crdate      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=10)     :: crtime      ! needed by F90 DATE_AND_TIME intrinsic
    character(len=5)      :: crzone      ! needed by F90 DATE_AND_TIME intrinsic
    integer, dimension(8) :: values      ! needed by F90 DATE_AND_TIME intrinsic

    logical :: has_std_latlon, has_ustag_latlon, has_vstag_latlon

  call error_handler(E_ERR,'nc_write_model_atts','routine not written',source,revision,revdate)

    if ( .not. module_initialized ) call static_init_model

    ierr = -1 ! assume things go poorly

    has_std_latlon   = .FALSE.
    has_ustag_latlon = .FALSE.
    has_vstag_latlon = .FALSE.

    if (allocated(cosmo_lonlat(1)%lon) .and. allocated(cosmo_lonlat(1)%lat)) has_std_latlon   = .TRUE.
    if (allocated(cosmo_lonlat(2)%lon) .and. allocated(cosmo_lonlat(2)%lat)) has_ustag_latlon = .TRUE.
    if (allocated(cosmo_lonlat(3)%lon) .and. allocated(cosmo_lonlat(3)%lat)) has_vstag_latlon = .TRUE.

    write(filename,*) 'ncFileID', ncFileID

    !-------------------------------------------------------------------------------
    ! make sure ncFileID refers to an open netCDF file,
    ! and then put into define mode.
    !-------------------------------------------------------------------------------

    call nc_check(nf90_Inquire(ncFileID,nDimensions,nVariables,nAttributes,unlimitedDimID),&
                                       'nc_write_model_atts', 'inquire '//trim(filename))
    call nc_check(nf90_Redef(ncFileID),'nc_write_model_atts',   'redef '//trim(filename))

    !-------------------------------------------------------------------------------
    ! We need the dimension ID for the number of copies/ensemble members, and
    ! we might as well check to make sure that Time is the Unlimited dimension.
    ! Our job is create the 'model size' dimension.
    !-------------------------------------------------------------------------------

    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='NMLlinelen', dimid=LineLenDimID), &
     'nc_write_model_atts','inq_dimid NMLlinelen')
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='copy', dimid=MemberDimID), &
     'nc_write_model_atts', 'copy dimid '//trim(filename))
    call nc_check(nf90_inq_dimid(ncid=ncFileID, name='time', dimid=  TimeDimID), &
     'nc_write_model_atts', 'time dimid '//trim(filename))

    if ( TimeDimID /= unlimitedDimId ) then
      write(string1,*)'Time Dimension ID ',TimeDimID, &
       ' should equal Unlimited Dimension ID',unlimitedDimID
!      call error_handler(E_ERR,'nc_write_model_atts', string1, source, revision, revdate)
    endif

    !-------------------------------------------------------------------------------
    ! Define the model size / state variable dimension / whatever ...
    !-------------------------------------------------------------------------------
    call nc_check(nf90_def_dim(ncid=ncFileID, name='StateVariable', len=model_size, &
     dimid = StateVarDimID),'nc_write_model_atts', 'state def_dim '//trim(filename))

    !-------------------------------------------------------------------------------
    ! Write Global Attributes
    !-------------------------------------------------------------------------------

     call DATE_AND_TIME(crdate,crtime,crzone,values)
     write(string1,'(''YYYY MM DD HH MM SS = '',i4,5(1x,i2.2))') &
      values(1), values(2), values(3), values(5), values(6), values(7)

     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'creation_date' ,string1 ), &
                   'nc_write_model_atts', 'creation put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_source'  ,source  ), &
                   'nc_write_model_atts', 'source put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revision',revision), &
                   'nc_write_model_atts', 'revision put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model_revdate' ,revdate ), &
                   'nc_write_model_atts', 'revdate put '//trim(filename))
     call nc_check(nf90_put_att(ncFileID, NF90_GLOBAL, 'model',  'cosmo' ), &
                   'nc_write_model_atts', 'model put '//trim(filename))

    !-------------------------------------------------------------------------------
    ! Here is the extensible part. The simplest scenario is to output the state vector,
    ! parsing the state vector into model-specific parts is complicated, and you need
    ! to know the geometry, the output variables (PS,U,V,T,Q,...) etc. We're skipping
    ! complicated part.
    !-------------------------------------------------------------------------------

    if ( output_1D_state_vector ) then

      !----------------------------------------------------------------------------
      ! Create a variable for the state vector
      !----------------------------------------------------------------------------

      ! Define the state vector coordinate variable and some attributes.
      call nc_check(nf90_def_var(ncid=ncFileID,name='StateVariable', xtype=nf90_int, &
                    dimids=StateVarDimID, varid=StateVarVarID), 'nc_write_model_atts', &
                    'statevariable def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarVarID,'long_name','State Variable ID'),&
                    'nc_write_model_atts','statevariable long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID, StateVarVarID, 'units','indexical'), &
                    'nc_write_model_atts', 'statevariable units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarVarID,'valid_range',(/ 1,model_size /)),&
                    'nc_write_model_atts', 'statevariable valid_range '//trim(filename))

      ! Define the actual (3D) state vector, which gets filled as time goes on ...
      call nc_check(nf90_def_var(ncid=ncFileID, name='state', xtype=nf90_real, &
                    dimids=(/StateVarDimID,MemberDimID,unlimitedDimID/),varid=StateVarID),&
                    'nc_write_model_atts','state def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,StateVarID,'long_name','model state or fcopy'),&
                    'nc_write_model_atts', 'state long_name '//trim(filename))

      ! Leave define mode so we can fill the coordinate variable.
      call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','state enddef '//trim(filename))

      ! Fill the state variable coordinate variable
      call nc_check(nf90_put_var(ncFileID, StateVarVarID, (/ (i,i=1,model_size) /) ), &
                    'nc_write_model_atts', 'state put_var '//trim(filename))

    else

      !----------------------------------------------------------------------------
      ! We need to output the prognostic variables.
      !----------------------------------------------------------------------------
      ! Define the new dimensions IDs
      !----------------------------------------------------------------------------

      findnxny : do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
          nx=state_vector_vars(ikind)%nx
          ny=state_vector_vars(ikind)%ny
          exit findnxny
        endif
      enddo findnxny

      findnz : do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then
          if ((state_vector_vars(ikind)%nz>1) .and. (ikind .ne. KIND_VERTICAL_VELOCITY)) then
            nz=state_vector_vars(ikind)%nz
            exit findnz
          endif
        endif
      enddo findnz

      call nc_check(nf90_def_dim(ncid=ncFileID, name='lon', len=nx, dimid = lonDimID), &
                     'nc_write_model_atts', 'lon def_dim '//trim(filename))
      call nc_check(nf90_def_dim(ncid=ncFileID, name='lat', len=ny, dimid = latDimID), &
                     'nc_write_model_atts', 'lat def_dim '//trim(filename))
      call nc_check(nf90_def_dim(ncid=ncFileID, name='lev', len=nz, dimid = levDimID), &
                     'nc_write_model_atts', 'lev def_dim '//trim(filename))
      call nc_check(nf90_def_dim(ncid=ncFileID, name='wlev', len=nz+1, dimid = wlevDimID), &
                     'nc_write_model_atts', 'lev def_dim '//trim(filename))

      if ( has_std_latlon ) then
      ! Standard Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='LON', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=lonVarID),&
                    'nc_write_model_atts', 'LON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'long_name', 'longitudes of grid'), &
                    'nc_write_model_atts', 'LON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'LON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'LON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  lonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'LON valid_range '//trim(filename))
      ! Standard Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='LAT', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=latVarID),&
                    'nc_write_model_atts', 'LAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'long_name', 'latitudes of grid'), &
                    'nc_write_model_atts', 'LAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'LAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'LAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  latVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'LAT valid_range '//trim(filename))
      endif


      if ( has_ustag_latlon ) then
      ! U Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='ULON', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=ulonVarID),&
                    'nc_write_model_atts', 'ULON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'long_name', 'longitudes for U-wind'), &
                    'nc_write_model_atts', 'ULON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'ULON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'ULON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'ULON valid_range '//trim(filename))
      ! U Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='ULAT', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=ulatVarID),&
                    'nc_write_model_atts', 'ULAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'long_name', 'latitudes for U-wind'), &
                    'nc_write_model_atts', 'ULAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'ULAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'ULAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  ulatVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'ULAT valid_range '//trim(filename))
      endif


      if ( has_vstag_latlon ) then
      ! V Grid Longitudes
      call nc_check(nf90_def_var(ncFileID,name='VLON', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=vlonVarID),&
                    'nc_write_model_atts', 'VLON def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'long_name', 'longitudes for V-wind'), &
                    'nc_write_model_atts', 'VLON long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'cartesian_axis', 'X'),  &
                    'nc_write_model_atts', 'VLON cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'VLON units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlonVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'VLON valid_range '//trim(filename))
      ! V Grid Latitudes
      call nc_check(nf90_def_var(ncFileID,name='VLAT', xtype=nf90_real, &
                    dimids=(/ lonDimID, latDimID /), varid=vlatVarID),&
                    'nc_write_model_atts', 'VLAT def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'long_name', 'latitudes for V-wind'), &
                    'nc_write_model_atts', 'VLAT long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'cartesian_axis', 'Y'),  &
                    'nc_write_model_atts', 'VLAT cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'units', 'degrees_east'), &
                    'nc_write_model_atts', 'VLAT units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  vlatVarID, 'valid_range', (/ -180.0_r8, 360.0_r8 /)), &
                    'nc_write_model_atts', 'VLAT valid_range '//trim(filename))
      endif

      ! Standard Z Levels
      call nc_check(nf90_def_var(ncFileID,name='LEV', xtype=nf90_real, &
                    dimids=(/ levDimID /), varid=levVarID),&
                    'nc_write_model_atts', 'LEV def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'long_name', 'standard hybrid model levels'), &
                    'nc_write_model_atts', 'LEV long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'cartesian_axis', 'Z'),  &
                    'nc_write_model_atts', 'LEV cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'units', 'model level'), &
                    'nc_write_model_atts', 'LEV units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  levVarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
                    'nc_write_model_atts', 'LEV valid_range '//trim(filename))

      ! W-wind Z Levels
      call nc_check(nf90_def_var(ncFileID,name='WLEV', xtype=nf90_real, &
                    dimids=(/ wlevDimID /), varid=wlevVarID),&
                    'nc_write_model_atts', 'WLEV def_var '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'long_name', 'standard model levels for W-wind'), &
                    'nc_write_model_atts', 'WLEV long_name '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'cartesian_axis', 'Z'),  &
                    'nc_write_model_atts', 'WLEV cartesian_axis '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'units', 'model level'), &
                    'nc_write_model_atts', 'WLEV units '//trim(filename))
      call nc_check(nf90_put_att(ncFileID,  wlevVarID, 'valid_range', (/ 1._r8,float(nz)+1._r8 /)), &
                    'nc_write_model_atts', 'WLEV valid_range '//trim(filename))

      ! DEBUG block to check shape of netCDF variables
      if ( debug > 0 .and. do_output() ) then
         write(*,*)'lon   dimid is ',lonDimID
         write(*,*)'lat   dimid is ',latDimID
         write(*,*)'lev   dimid is ',levDimID
         write(*,*)'wlev  dimid is ',wlevDimID
         write(*,*)'unlim dimid is ',unlimitedDimID
         write(*,*)'copy  dimid is ',MemberDimID
      endif

      do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then

          string1 = trim(filename)//' '//trim(state_vector_vars(ikind)%varname_short)

          dims(1)=lonDimID
          dims(2)=latDimID

          idim=3
          if (state_vector_vars(ikind)%nz>1) then
            dims(idim)=levDimID
            if (ikind==KIND_VERTICAL_VELOCITY) then
              wlevs(1:nz+1)=state_vector_vars(ikind)%vertical_level(1:nz+1)
              dims(idim)=wlevDimID
            else
              levs(1:nz)=state_vector_vars(ikind)%vertical_level(1:nz)
            endif
            idim=idim+1
          endif

          ! Create a dimension for the ensemble
          dims(idim) = memberDimID
          idim=idim+1

          ! Put ensemble member dimension here
          dims(idim) = unlimitedDimID
          ndims=idim

          ! check shape of netCDF variables
          if ( debug > 0 .and. do_output() ) &
          write(*,*)trim(state_vector_vars(ikind)%varname_short),' has netCDF dimIDs ',dims(1:ndims)

          call nc_check(nf90_def_var(ncid=ncFileID, name=trim(state_vector_vars(ikind)%varname_short), xtype=nf90_real, &
                        dimids = dims(1:ndims), varid=VarID),&
                        'nc_write_model_atts', trim(string1)//' def_var' )

          call nc_check(nf90_put_att(ncFileID, VarID, 'long_name', trim(state_vector_vars(ikind)%varname_long)), &
                        'nc_write_model_atts', trim(string1)//' put_att long_name' )

          write(ckind,'(I6)') ikind
          call nc_check(nf90_put_att(ncFileID, VarID, 'DART_kind', trim(ckind)), &
                        'nc_write_model_atts', trim(string1)//' put_att dart_kind' )

          call nc_check(nf90_put_att(ncFileID, VarID, 'units', trim(state_vector_vars(ikind)%units)), &
                        'nc_write_model_atts', trim(string1)//' put_att units' )

        endif
      enddo

      ! Leave define mode so we can fill the coordinate variable.
      call nc_check(nf90_enddef(ncfileID),'nc_write_model_atts','prognostic enddef '//trim(filename))

      !----------------------------------------------------------------------------
      ! Fill the coordinate variables - reshape the 1D arrays to the 2D shape
      !----------------------------------------------------------------------------
      allocate(data2d(nx,ny))

      if (has_std_latlon) then
      data2d = reshape(cosmo_lonlat(1)%lon, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, lonVarID, data2d), &
                    'nc_write_model_atts', 'LON put_var '//trim(filename))

      data2d = reshape(cosmo_lonlat(1)%lat, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, latVarID, data2d ), &
                    'nc_write_model_atts', 'LAT put_var '//trim(filename))
      endif

      if (has_ustag_latlon) then
      data2d = reshape(cosmo_lonlat(2)%lon, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, ulonVarID, data2d ), &
                    'nc_write_model_atts', 'ULON put_var '//trim(filename))

      data2d = reshape(cosmo_lonlat(2)%lat, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, ulatVarID, data2d ), &
                    'nc_write_model_atts', 'ULAT put_var '//trim(filename))
      endif

      if (has_vstag_latlon) then
      data2d = reshape(cosmo_lonlat(3)%lon, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, vlonVarID, data2d ), &
                    'nc_write_model_atts', 'VLON put_var '//trim(filename))

      data2d = reshape(cosmo_lonlat(3)%lat, (/ nx, ny /) )
      call nc_check(nf90_put_var(ncFileID, vlatVarID, data2d ), &
                    'nc_write_model_atts', 'VLAT put_var '//trim(filename))
      endif

      deallocate(data2d)

      call nc_check(nf90_put_var(ncFileID, levVarID, levs(1:nz) ), &
                    'nc_write_model_atts', 'LEV put_var '//trim(filename))

      call nc_check(nf90_put_var(ncFileID, wlevVarID, wlevs(1:nz+1) ), &
                    'nc_write_model_atts', 'WLEV put_var '//trim(filename))

    endif

    !-------------------------------------------------------------------------------
    ! Flush the buffer and leave netCDF file open
    !-------------------------------------------------------------------------------
    call nc_check(nf90_sync(ncFileID), 'nc_write_model_atts', 'atts sync')

    ierr = 0 ! If we got here, things went well.

  end function nc_write_model_atts


!------------------------------------------------------------------------
!>


  function nc_write_model_vars( ncFileID, state_vec, copyindex, timeindex ) result (ierr)
    !------------------------------------------------------------------
    ! TJH 24 Oct 2006 -- Writes the model variables to a netCDF file.
    !
    ! TJH 29 Jul 2003 -- for the moment, all errors are fatal, so the
    ! return code is always '0 == normal', since the fatal errors stop execution.
    !
    ! For the lorenz_96 model, each state variable is at a separate location.
    ! that's all the model-specific attributes I can think of ...
    !
    ! assim_model_mod:init_diag_output uses information from the location_mod
    !     to define the location dimension and variable ID. All we need to do
    !     is query, verify, and fill ...
    !
    ! Typical sequence for adding new dimensions,variables,attributes:
    ! NF90_OPEN             ! open existing netCDF dataset
    !    NF90_redef         ! put into define mode
    !    NF90_def_dim       ! define additional dimensions (if any)
    !    NF90_def_var       ! define variables: from name, type, and dims
    !    NF90_put_att       ! assign attribute values
    ! NF90_ENDDEF           ! end definitions: leave define mode
    !    NF90_put_var       ! provide values for variable
    ! NF90_CLOSE            ! close: save updated netCDF dataset

    integer,                intent(in) :: ncFileID      ! netCDF file identifier
    real(r8), dimension(:), intent(in) :: state_vec
    integer,                intent(in) :: copyindex
    integer,                intent(in) :: timeindex
    integer                            :: ierr          ! return value of function

    integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs, mystart, mycount
    character(len=NF90_MAX_NAME)          :: varname
    integer :: i,ikind, VarID, ncNdims, dimlen,ndims,vardims(3)
    integer :: TimeDimID, CopyDimID

    real(r8), allocatable, dimension(:,:)     :: data_2d_array
    real(r8), allocatable, dimension(:,:,:)   :: data_3d_array

    character(len=128) :: filename

  call error_handler(E_ERR,'nc_write_model_vars','routine not written',source,revision,revdate)

    if ( .not. module_initialized ) call static_init_model

    ierr = -1 ! assume things go poorly

    !--------------------------------------------------------------------
    ! we only have a netcdf handle here so we do not know the filename
    ! or the fortran unit number.  but construct a string with at least
    ! the netcdf handle, so in case of error we can trace back to see
    ! which netcdf file is involved.
    !--------------------------------------------------------------------

    write(filename,*) 'ncFileID', ncFileID

    !-------------------------------------------------------------------------------
    ! make sure ncFileID refers to an open netCDF file,
    !-------------------------------------------------------------------------------

    call nc_check(nf90_inq_dimid(ncFileID, 'copy', dimid=CopyDimID), &
     'nc_write_model_vars', 'inq_dimid copy '//trim(filename))

    call nc_check(nf90_inq_dimid(ncFileID, 'time', dimid=TimeDimID), &
     'nc_write_model_vars', 'inq_dimid time '//trim(filename))

    if ( output_1D_state_vector ) then

      call nc_check(NF90_inq_varid(ncFileID, 'state', VarID), &
       'nc_write_model_vars', 'state inq_varid '//trim(filename))
      call nc_check(NF90_put_var(ncFileID,VarID,state_vec,start=(/1,copyindex,timeindex/)),&
       'nc_write_model_vars', 'state put_var '//trim(filename))

    else

      !----------------------------------------------------------------------------
      ! We need to process the prognostic variables.
      !----------------------------------------------------------------------------

      do ikind=1,n_max_kinds
        if (state_vector_vars(ikind)%is_present) then

          varname = trim(state_vector_vars(ikind)%varname_short)
          string1 = trim(filename)//' '//trim(varname)

          ! Ensure netCDF variable is conformable with progvar quantity.
          ! The TIME and Copy dimensions are intentionally not queried
          ! by looping over the dimensions stored in the progvar type.

          call nc_check(nf90_inq_varid(ncFileID, varname, VarID), &
                        'nc_write_model_vars', 'inq_varid '//trim(string1))

          call nc_check(nf90_inquire_variable(ncFileID,VarId,dimids=dimIDs,ndims=ncNdims), &
                        'nc_write_model_vars', 'inquire '//trim(string1))

          mystart(:)=1
          mycount(:)=1

          if (state_vector_vars(ikind)%nz==1) then
            ndims=2
          else
            ndims=3
          endif

          vardims(1)=state_vector_vars(ikind)%nx
          vardims(2)=state_vector_vars(ikind)%ny
          vardims(3)=state_vector_vars(ikind)%nz

          DimCheck : do i = 1,ndims

            write(string1,'(a,i2,A)') 'inquire dimension ',i,trim(varname)
            call nc_check(nf90_inquire_dimension(ncFileID, dimIDs(i), len=dimlen), &
             'nc_write_model_vars', string1)

            if ( dimlen /= vardims(i) ) then
              write(string1,*) trim(varname),' dim/dimlen ',i,dimlen,' not ',vardims(i)
              write(string2,*)' but it should be.'
              call error_handler(E_ERR, 'nc_write_model_vars', string1, &
                              source, revision, revdate, text2=string2)
            endif

            mycount(i) = dimlen

          enddo DimCheck

          where(dimIDs == CopyDimID) mystart = copyindex
          where(dimIDs == CopyDimID) mycount = 1
          where(dimIDs == TimeDimID) mystart = timeindex
          where(dimIDs == TimeDimID) mycount = 1

          if (ndims==2) then
            allocate(data_2d_array(vardims(1),vardims(2)))
            call sv_to_field(data_2d_array,state_vec,state_vector_vars(ikind))
            call nc_check(nf90_put_var(ncFileID, VarID, data_2d_array, &
                          start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                          'nc_write_model_vars', 'put_var '//trim(string2))
            deallocate(data_2d_array)

          elseif (ndims==3) then
            allocate(data_3d_array(vardims(1),vardims(2),vardims(3)))
            call sv_to_field(data_3d_array,state_vec,state_vector_vars(ikind))
            call nc_check(nf90_put_var(ncFileID, VarID, data_3d_array, &
                          start = mystart(1:ncNdims), count=mycount(1:ncNdims)), &
                          'nc_write_model_vars', 'put_var '//trim(string2))
            deallocate(data_3d_array)

          else
             write(string1, *) 'no support for data array of dimension ', ncNdims
             call error_handler(E_ERR,'nc_write_model_vars', string1, &
                           source,revision,revdate)
          endif

        endif

      enddo

    endif

    return

  end function nc_write_model_vars


!------------------------------------------------------------------------
!>


  subroutine pert_model_state(state, pert_state, interf_provided)
    !------------------------------------------------------------------
    ! Perturbs a model state for generating initial ensembles.
    ! The perturbed state is returned in pert_state.
    ! A model may choose to provide a NULL INTERFACE by returning
    ! .false. for the interf_provided argument. This indicates to
    ! the filter that if it needs to generate perturbed states, it
    ! may do so by adding a perturbation to each model state
    ! variable independently. The interf_provided argument
    ! should be returned as .true. if the model wants to do its own
    ! perturbing of states.
    !------------------------------------------------------------------
    ! Currently only implemented as rondom perturbations
    !------------------------------------------------------------------

    real(r8), intent(in)  :: state(:)
    real(r8), intent(out) :: pert_state(:)
    logical,  intent(out) :: interf_provided

    real(r8)              :: stddev,mean

    integer               :: ikind,ilevel,i,istart,iend
    logical, save         :: random_seq_init = .false.

  call error_handler(E_ERR,'pert_model_state','routine not written',source,revision,revdate)
    if ( .not. module_initialized ) call static_init_model

    interf_provided = .true.

    ! Initialize my random number sequence (no seed is submitted here!)
    if(.not. random_seq_init) then
      call init_random_seq(random_seq)
      random_seq_init = .true.
    endif

    ! add some uncertainty to every state vector element
    do ikind=1,size(state_vector_vars)
      if (state_vector_vars(ikind)%is_present) then
        do ilevel=1,state_vector_vars(ikind)%nz
          istart=state_vector_vars(ikind)%state_vector_sindex(ilevel)
          iend=istart+(state_vector_vars(ikind)%nx*state_vector_vars(ikind)%ny)-1

          mean=sum(abs(state(istart:iend)))/float(iend-istart+1)
          stddev=sqrt(sum((state(istart:iend)-mean)**2))/float(iend-istart+1)

          do i=istart,iend
            pert_state(i) = random_gaussian(random_seq, state(i),model_perturbation_amplitude*stddev)
          enddo
          if ((ikind==KIND_SPECIFIC_HUMIDITY) .or. &
              (ikind==KIND_CLOUD_LIQUID_WATER) .or. &
              (ikind==KIND_CLOUD_ICE)) then
            where (pert_state(istart:iend)<0.)
              pert_state(istart:iend)=0.
            end where
          endif
        enddo
      endif
    enddo

    return

  end subroutine pert_model_state


!------------------------------------------------------------------------
!>


  subroutine ens_mean_for_model(filter_ens_mean)

    real(r8), dimension(:), intent(in) :: filter_ens_mean

  call error_handler(E_ERR,'ens_mean_for_model','routine not written',source,revision,revdate)

    if ( .not. module_initialized ) call static_init_model

    allocate(ens_mean(1:model_size))
    ens_mean(:) = filter_ens_mean(:)

!  write(string1,*) 'COSMO has no ensemble mean in storage.'
!  call error_handler(E_ERR,'ens_mean_for_model',string1,source,revision,revdate)

  end subroutine ens_mean_for_model


!------------------------------------------------------------------------
!>


subroutine set_allowed_state_vector_vars()
  ! set the information on which variables should go into the state vector

  is_allowed_state_vector_var(:)=.FALSE.
  is_allowed_non_state_var(:)=.FALSE.

  allowed_state_vector_vars(1)=KIND_U_WIND_COMPONENT
   is_allowed_state_vector_var(KIND_U_WIND_COMPONENT)=.TRUE.

  allowed_state_vector_vars(2)=KIND_V_WIND_COMPONENT
   is_allowed_state_vector_var(KIND_V_WIND_COMPONENT)=.TRUE.

  allowed_state_vector_vars(3)=KIND_VERTICAL_VELOCITY
   is_allowed_state_vector_var(KIND_VERTICAL_VELOCITY)=.TRUE.

  allowed_state_vector_vars(4)=KIND_TEMPERATURE
   is_allowed_state_vector_var(KIND_TEMPERATURE)=.TRUE.

  allowed_state_vector_vars(5)=KIND_PRESSURE
   is_allowed_state_vector_var(KIND_PRESSURE)=.TRUE.

  allowed_state_vector_vars(6)=KIND_SPECIFIC_HUMIDITY
   is_allowed_state_vector_var(KIND_SPECIFIC_HUMIDITY)=.TRUE.

  allowed_state_vector_vars(7)=KIND_CLOUD_LIQUID_WATER
   is_allowed_state_vector_var(KIND_CLOUD_LIQUID_WATER)=.TRUE.

  allowed_state_vector_vars(8)=KIND_CLOUD_ICE
   is_allowed_state_vector_var(KIND_CLOUD_ICE)=.TRUE.

  ! set the information which variables are needed but will not go into the state vector
  allowed_non_state_vars(1)=KIND_SURFACE_ELEVATION
   is_allowed_non_state_var(KIND_SURFACE_ELEVATION)=.TRUE.
  allowed_non_state_vars(2)=KIND_SURFACE_GEOPOTENTIAL
   is_allowed_non_state_var(KIND_SURFACE_GEOPOTENTIAL)=.TRUE.
  allowed_non_state_vars(3)=KIND_PRESSURE_PERTURBATION
   is_allowed_non_state_var(KIND_PRESSURE_PERTURBATION)=.TRUE.

  return

end subroutine set_allowed_state_vector_vars


!------------------------------------------------------------------------
!>


  function ll_to_xyz_vector(lon,lat) RESULT (xyz)

    ! Passed variables

    real(r8),allocatable :: xyz(:,:)      ! result: x,z,y-coordinates
    real(r8),intent(in)  :: lat(:),lon(:) ! input:  lat/lon coordinates in degrees

    real(r8)             :: radius
    integer              :: n

    ! define output vector size to be the same as the input vector size
    ! second dimension (3) is x,y,z

    n=SIZE(lat,1)
    ALLOCATE(xyz(1:n,1:3))

    ! as we are interested in relative distances we set the radius to 1 - may be changed later

    radius=1.0_r8

    ! caclulate the x,y,z-coordinates

    xyz(1:n,1)=radius*sin(lat(1:n)*deg2rad)*cos(lon(1:n)*deg2rad)
    xyz(1:n,2)=radius*sin(lat(1:n)*deg2rad)*sin(lon(1:n)*deg2rad)
    xyz(1:n,3)=radius*cos(lat(1:n)*deg2rad)

    return
  end function ll_to_xyz_vector


!------------------------------------------------------------------------
!>


  function ll_to_xyz_single(lon,lat) result (xyz)

    ! Passed variables

    real(r8)             :: xyz(1:3) ! result: x,z,y-coordinates
    real(r8),intent(in)  :: lat,lon  ! input:  lat/lon coordinates in degrees

    real(r8)             :: radius

    ! as we are interested in relative distances we set the radius to 1 - may be changed later

    radius=1.0_r8

    ! caclulate the x,y,z-coordinates

    xyz(1)=radius*sin(lat*deg2rad)*cos(lon*deg2rad)
    xyz(2)=radius*sin(lat*deg2rad)*sin(lon*deg2rad)
    xyz(3)=radius*cos(lat*deg2rad)

    return
  end function ll_to_xyz_single


!------------------------------------------------------------------------
!>


  subroutine get_enclosing_grid_box(p,g,n,nx,ny,b,bw)

    integer,intent(in)   :: n,nx,ny
    real(r8),intent(in)  :: p(1:3),g(1:n,1:3)
    integer,intent(out)  :: b(1:2,1:2)
    real(r8),intent(out) :: bw(1:2,1:2)

!    real(r8)             :: work(1:nx,1:ny,1:3),dist(1:nx,1:ny),boxdist(1:2,1:2)
    real(r8)             :: work(1:nx+2,1:ny+2,1:3),dist(1:nx+2,1:ny+2),boxdist(1:2,1:2)
    integer              :: i,j,minidx(2),boxidx(2),xb,yb

    real(r8) :: sqrt2

    sqrt2 = sqrt(2.0_r8)

    work(2:nx+1,2:ny+1,1:3)=RESHAPE( g, (/ nx,ny,3 /))

    do i=2,nx+1
      work(i,   1,1:3)=work(i,   2,1:3)-(work(i, 3,1:3)-work(i,   2,1:3))
      work(i,ny+2,1:3)=work(i,ny+1,1:3)-(work(i,ny,1:3)-work(i,ny+1,1:3))
    enddo

    do j=2,ny+1
      work(   1,j,1:3)=work(   2,j,1:3)-(work( 3,j,1:3)-work(   2,j,1:3))
      work(nx+2,j,1:3)=work(nx+1,j,1:3)-(work(nx,j,1:3)-work(nx+1,j,1:3))
    enddo

    work(   1,   1,1:3) = work(   2,   2,1:3) - 0.5_r8*(sqrt2*(work(   2,   2,1:3)-work(   1,   2,1:3)) + sqrt2*(work(   2,   2,1:3)-work(   2,   1,1:3)))
    work(   1,ny+2,1:3) = work(   2,ny+1,1:3) - 0.5_r8*(sqrt2*(work(   2,ny+1,1:3)-work(   1,ny+1,1:3)) + sqrt2*(work(   2,ny+1,1:3)-work(   2,ny+2,1:3)))
    work(nx+2,   1,1:3) = work(nx+1,   2,1:3) - 0.5_r8*(sqrt2*(work(nx+1,   2,1:3)-work(nx+2,   2,1:3)) + sqrt2*(work(nx+1,   2,1:3)-work(nx+1,   1,1:3)))
    work(nx+2,ny+2,1:3) = work(nx+1,ny+1,1:3) - 0.5_r8*(sqrt2*(work(nx+1,ny+1,1:3)-work(nx+2,ny+1,1:3)) + sqrt2*(work(nx+1,ny+1,1:3)-work(nx+1,ny+2,1:3)))

    do i=1,nx+2
    do j=1,ny+2
        dist(i,j)=sqrt(sum((work(i,j,:)-p(:))**2))
    enddo
    enddo

    minidx(:)=minloc(dist)

    ! watch for out of area values

    if (minidx(1)==1 .or. minidx(1)==(nx+2) .or. minidx(2)==1 .or. minidx(2)==(ny+2)) then
      b(:,:)=-1
      return
    endif


    do i=0,1
    do j=0,1
        boxdist(i+1,j+1)=sum(dist(minidx(1)+i-1:minidx(1)+i,minidx(2)+j-1:minidx(2)+j))
    enddo
    enddo

    boxidx=minloc(boxdist)-1

    xb=minidx(1)+(2*(boxidx(1)-0.5))
    yb=minidx(2)+(2*(boxidx(2)-0.5))

    if (xb==1 .or. xb==(nx+2) .or. yb==1 .or. yb==(ny+2)) then
      b(:,:)=-1
      return
    else
      do i=1,2
      do j=1,2
          b(i,j)=((minidx(2)+(j-1)*(boxidx(2)-0.5)*2)*ny)+(minidx(1)+(i-1)*(2*(boxidx(1)-0.5)))
      enddo
      enddo

      do i=1,2
      do j=1,2
          boxdist(i,j)=dist(mod(b(i,j),ny),b(i,j)/ny)
      enddo
      enddo

      bw(:,:)=1./boxdist(:,:)
!      bw(:,:)=(((1.-boxdist(:,:))/(1.1*maxval(boxdist)))**2)/((boxdist(:,:)/(1.1*maxval(boxdist)))**2)
      bw=bw/sum(bw)
      b(:,:)=b(:,:)-1
    endif

  end subroutine get_enclosing_grid_box


!------------------------------------------------------------------------
!>


  subroutine get_enclosing_grid_box_lonlat(lon,lat,p,n,nx,ny,b,bw)

    integer, intent(in)  :: n,nx,ny
    real(r8),intent(in)  :: p(1:2),lon(1:n),lat(1:n)
    integer, intent(out) :: b(1:2,1:2)
    real(r8),intent(out) :: bw(1:2,1:2)

!    real(r8)            :: work(1:nx,1:ny,1:3),dist(1:nx,1:ny),boxdist(1:2,1:2)
    real(r8)             :: work(1:nx+2,1:ny+2,1:2),dist(1:nx+2,1:ny+2),boxdist(1:2,1:2),pw(2)

    integer  :: i,j,minidx(2),boxidx(2),xb,yb,bx(2,2),by(2,2)
    real(r8) :: sqrt2

    sqrt2 = sqrt(2.0_r8)

    work(2:nx+1,2:ny+1,1)=reshape(lon,(/ nx,ny /))*deg2rad
    work(2:nx+1,2:ny+1,2)=reshape(lat,(/ nx,ny /))*deg2rad
    pw=p*deg2rad

    do i=2,nx+1
      work(i,   1,1:2)=work(i,   2,1:2)-(work(i, 3,1:2)-work(i,   2,1:2))
      work(i,ny+2,1:2)=work(i,ny+1,1:2)-(work(i,ny,1:2)-work(i,ny+1,1:2))
    enddo

    do j=2,ny+1
      work(   1,j,1:2)=work(   2,j,1:2)-(work( 3,j,1:2)-work(   2,j,1:2))
      work(nx+2,j,1:2)=work(nx+1,j,1:2)-(work(nx,j,1:2)-work(nx+1,j,1:2))
    enddo

    work(   1,   1,1:2) = work(   2,   2,1:2) - 0.5_r8*(sqrt2*(work(   2,   2,1:2)-work(   1,   2,1:2))+sqrt2*(work(   2,   2,1:2)-work(   2,   1,1:2)))
    work(   1,ny+2,1:2) = work(   2,ny+1,1:2) - 0.5_r8*(sqrt2*(work(   2,ny+1,1:2)-work(   1,ny+1,1:2))+sqrt2*(work(   2,ny+1,1:2)-work(   2,ny+2,1:2)))
    work(nx+2,   1,1:2) = work(nx+1,   2,1:2) - 0.5_r8*(sqrt2*(work(nx+1,   2,1:2)-work(nx+2,   2,1:2))+sqrt2*(work(nx+1,   2,1:2)-work(nx+1,   1,1:2)))
    work(nx+2,ny+2,1:2) = work(nx+1,ny+1,1:2) - 0.5_r8*(sqrt2*(work(nx+1,ny+1,1:2)-work(nx+2,ny+1,1:2))+sqrt2*(work(nx+1,ny+1,1:2)-work(nx+1,ny+2,1:2)))

    do i=1,nx+2
    do j=1,ny+2
!      dist(i,j)=sqrt(sum((work(i,j,:)-p(:))**2))
       dist(i,j) = 6173.0_r8*acos(cos(work(i,j,2)-pw(2))-cos(work(i,j,2))*cos(pw(2))*(1-cos(work(i,j,1)-pw(1))))
    enddo
    enddo

    minidx(:)=minloc(dist)

    ! watch for out of area values

    if (minidx(1)==1 .or. minidx(1)==(nx+2) .or. minidx(2)==1 .or. minidx(2)==(ny+2)) then
      b(:,:)=-1
      return
    endif

!   open(21,file='/daten02/jkeller/testbox.bin',form='unformatted')
!   iunit = open_file('testbox.bin',form='unformatted',action='write')
!   write(iunit) nx
!   write(iunit) ny

    do i=0,1
    do j=0,1
        boxdist(i+1,j+1)=sum(dist(minidx(1)+i-1:minidx(1)+i,minidx(2)+j-1:minidx(2)+j))/4.0_r8
!       write(*,'(4(I5))') minidx(1)+i-1,minidx(1)+i,minidx(2)+j-1,minidx(2)+j
!       write(iunit) (minidx(2)+j-1),minidx(1)+i-1,&
!                 (minidx(2)+j-1),minidx(1)+i,&
!                 (minidx(2)+j),minidx(1)+i-1,&
!                 (minidx(2)+j),minidx(1)+i
!       write(iunit) boxdist(i+1,j+1)
    enddo
    enddo

    boxidx=minloc(boxdist)-1

    xb=minidx(1)+(2*(boxidx(1)-0.5_r8))
    yb=minidx(2)+(2*(boxidx(2)-0.5_r8))

    if (xb==1 .or. xb==(nx+2) .or. yb==1 .or. yb==(ny+2)) then
      b(:,:)=-1
      return
    else
      do i=1,2
      do j=1,2
          bx(i,j)=minidx(1)+(i-1)*(2*(boxidx(1)-0.5_r8))
          by(i,j)=minidx(2)+(j-1)*(2*(boxidx(2)-0.5_r8))
      enddo
      enddo

      do i=1,2
      do j=1,2
          boxdist(i,j)=dist(bx(i,j),by(i,j))
      enddo
      enddo

      bw(:,:)=1.0_r8/boxdist(:,:)
      bw=bw/sum(bw)
      bx=bx-1
      by=by-1
      b(:,:)=(by-1)*nx+bx
    endif

    return

  end subroutine get_enclosing_grid_box_lonlat


!------------------------------------------------------------------------
!>


  subroutine bilinear_interpolation(bv,blo,bla,p,v)

    ! Passed variables

    real(r8),intent(in)  :: bv(2,2),blo(2,2),bla(2,2)
    real(r8),intent(in)  :: p(3)
    real(r8),intent(out) :: v

    ! Local storage

    real(r8)             :: x1,lo1,la1
    real(r8)             :: x2,lo2,la2
    real(r8)             :: d1,d2,d

!    write(*,'(3(F8.5,1X))') bv(1,1),blo(1,1),bla(1,1)
!    write(*,'(3(F8.5,1X))') bv(2,1),blo(2,1),bla(2,1)

    call linear_interpolation(p(1),p(2),bv(1,1),blo(1,1),bla(1,1),&
                                        bv(2,1),blo(2,1),bla(2,1),&
                                        x1,lo1,la1)

!    write(*,'(3(F8.5,1X))') x1,lo1,la1

    call linear_interpolation(p(1),p(2),bv(1,2),blo(1,2),bla(1,2),&
                                        bv(2,2),blo(2,2),bla(2,2),&
                                        x2,lo2,la2)

!    write(*,'(3(F8.5,1X))') x2,lo2,la2

    d1=sqrt((lo1-p(1))**2+(la1-p(2))**2)
    d2=sqrt((lo2-p(1))**2+(la2-p(2))**2)
    d =sqrt((lo1-lo2 )**2+(la1-la2 )**2)

    v=(1.0_r8-(d1/d))*x1+(1.0_r8-(d2/d))*x2

    return

  end subroutine bilinear_interpolation


!------------------------------------------------------------------------
!>


  subroutine linear_interpolation(lop,lap,x1,lo1,la1,x2,lo2,la2,x,lo,la)

    real(r8),intent(in)  :: lo1,lo2,la1,la2,x1,x2,lop,lap
    real(r8),intent(out) :: lo,la,x

    real(r8)             :: m1,m2,n1,n2,d1,d2,d,mylo1,mylo2,mylop,w1,w2

    mylo1=lo1
    mylo2=lo2
    mylop=lop

    if (lo1>180.0_r8) mylo1=lo1-360.0_r8
    if (lo2>180.0_r8) mylo2=lo2-360.0_r8
    if (lop>180.0_r8) mylop=lop-360.0_r8

    m1=(la2-la1)/(mylo2-mylo1)
    if (m1 .ne. 0.0_r8) then
      n1=la1-mylo1*m1
      m2=-1.0_r8/m1
      n2=lap-mylop*m2
      lo=(n2-n1)/(m1-m2)
      la=lo*m1+n1
      d1=sqrt((mylo1-lo)**2+(la1-la)**2)
      d2=sqrt((mylo2-lo)**2+(la2-la)**2)
      d =sqrt((mylo1-mylo2)**2+(la1-la2)**2)
    else
      la=la1
      lo=mylop

      d1=sqrt((mylo1-lo)**2+(la1-la)**2)
      d2=sqrt((mylo2-lo)**2+(la2-la)**2)
      d =sqrt((mylo1-mylo2)**2+(la1-la2)**2)
    endif

    if (lo < 0.0_r8) lo=lo+360.0_r8

    w1=abs(1.0_r8-(d1/d))
    w2=abs(1.0_r8-(d2/d))
    x=w1*x1+w2*x2

    return

  end subroutine linear_interpolation


!------------------------------------------------------------------------
!>


  subroutine get_vertical_boundaries(hb,hw,otype,vcs,p,b,w,istatus)

    real(r8), intent(in)  :: hw(2,2),p,vcs
    integer,  intent(in)  :: hb(2,2),otype
    integer,  intent(out) :: b(2),istatus
    real(r8), intent(out) :: w(2)

    integer               :: k,nlevel,x1,x2,x3,x4,y1,y2,y3,y4
    real(r8)              :: u,l
    real(r8),allocatable  :: klevel(:),hlevel(:),plevel(:)

    b(:)=-1

    ! coordinate system not implemented
    if ( (nint(vcs) == VERTISUNDEF)        .or. &
         (nint(vcs) == VERTISSURFACE)      .or. &
         (nint(vcs) == VERTISSCALEHEIGHT) ) then
      istatus=19
      return
    endif

! TJH    write(*,*)'non_state_data%pfl min max ',minval(non_state_data%pfl),maxval(non_state_data%pfl)
! TJH    write(*,*)' mean is ',sum(non_state_data%pfl)/(665.0_r8*657.0_r8*40.0_r8)

    x1 = mod(hb(1,1),size(non_state_data%pfl,1))
    x2 = mod(hb(2,1),size(non_state_data%pfl,1))
    x3 = mod(hb(1,2),size(non_state_data%pfl,1))
    x4 = mod(hb(2,2),size(non_state_data%pfl,1))
    y1 =     hb(1,1)/size(non_state_data%pfl,1)
    y2 =     hb(2,1)/size(non_state_data%pfl,1)
    y3 =     hb(1,2)/size(non_state_data%pfl,1)
    y4 =     hb(2,2)/size(non_state_data%pfl,1)

! TJH    write(*,*)'hb is ',hb
! TJH    write(*,*)'x  is ',x1,x2,x3,x4
! TJH    write(*,*)'y  is ',y1,y2,y3,y4
! TJH    write(*,*)'hw is ',hw

    if (otype .ne. KIND_VERTICAL_VELOCITY) then
      ! The variable exists on the 'full' levels
      nlevel=non_state_data%nfl
      allocate(klevel(1:nlevel))
      allocate(hlevel(1:nlevel))
      allocate(plevel(1:nlevel))
      klevel=state_vector_vars(otype)%vertical_level(:)
      hlevel=hw(1,1)*non_state_data%hfl(x1,y1,:)+&
             hw(2,1)*non_state_data%hfl(x2,y2,:)+&
             hw(1,2)*non_state_data%hfl(x3,y3,:)+&
             hw(2,2)*non_state_data%hfl(x4,y4,:)
      plevel=hw(1,1)*non_state_data%pfl(x1,y1,:)+&
             hw(2,1)*non_state_data%pfl(x2,y2,:)+&
             hw(1,2)*non_state_data%pfl(x3,y3,:)+&
             hw(2,2)*non_state_data%pfl(x4,y4,:)

! TJH write(*,*)non_state_data%pfl(x1,y1,:)
! TJH write(*,*)non_state_data%pfl(x2,y2,:)
! TJH write(*,*)non_state_data%pfl(x3,y3,:)
! TJH write(*,*)non_state_data%pfl(x4,y4,:)

    else
      ! The variable exists on the 'half' levels
      nlevel=non_state_data%nhl
      allocate(klevel(1:nlevel))
      allocate(hlevel(1:nlevel))
      allocate(plevel(1:nlevel))
      klevel=state_vector_vars(otype)%vertical_level(:)
      hlevel=hw(1,1)*non_state_data%hhl(x1,y1,:)+&
             hw(2,1)*non_state_data%hhl(x2,y2,:)+&
             hw(1,2)*non_state_data%hhl(x3,y3,:)+&
             hw(2,2)*non_state_data%hhl(x4,y4,:)
      plevel=hw(1,1)*non_state_data%phl(x1,y1,:)+&
             hw(2,1)*non_state_data%phl(x2,y2,:)+&
             hw(1,2)*non_state_data%phl(x3,y3,:)+&
             hw(2,2)*non_state_data%phl(x4,y4,:)
    endif

    u = -1.0_r8
    l = -1.0_r8

    do k=1,nlevel-1

      ! Find the bounding levels for the respective coordinate system
      if (nint(vcs) == VERTISLEVEL) then
        u=klevel(k+1)
        l=klevel(k)
      endif
      if (nint(vcs) == VERTISPRESSURE) then
      ! write(*,*)' vert is pressure '
      ! write(*,*)'plevel is ',plevel
        u=plevel(k+1)
        l=plevel(k)
      endif
      if (nint(vcs) == VERTISHEIGHT) then
      ! write(*,*)' vert is height '
      ! write(*,*)'hlevel is ',hlevel
        u=hlevel(k+1)
        l=hlevel(k)
      endif

! TJH write(*,*)'u p l',u,p,l

      if (u>=p .and. l<=p) then
        b(1)=k
        b(2)=k+1
        w(1)=1.0_r8-(p-l)/(u-l)
        w(2)=1.0_r8-(u-p)/(u-l)
        return
      endif

    enddo

    istatus=16 ! out of domain
    return
  end subroutine get_vertical_boundaries


!------------------------------------------------------------------------
!>


  subroutine sv_to_field_2d(f,x,v)

    real(r8),intent(out)                :: f(:,:)
    real(r8),intent(in)                 :: x(:)
    type(dart_variable_info),intent(in) :: v

    integer                             :: is,ie

    is=v%state_vector_sindex(1)
    ie=is+v%nx*v%ny-1
    f(:,:) = reshape(x(is:ie),(/ v%nx,v%ny /))

    return

  end subroutine sv_to_field_2d


!------------------------------------------------------------------------
!>


  subroutine sv_to_field_3d(f,x,v)

    real(r8),intent(out)                :: f(:,:,:)
    real(r8),intent(in)                 :: x(:)
    type(dart_variable_info),intent(in) :: v

    integer                             :: is,ie,iz

    do iz=1,v%nz
      is=v%state_vector_sindex(iz)
      ie=is+v%nx*v%ny-1
      f(:,:,iz) = reshape(x(is:ie),(/ v%nx,v%ny /))
    enddo

    return

  end subroutine sv_to_field_3d


!------------------------------------------------------------------------
!>


  function get_state_time() result (time)
    type(time_type) :: time

  call error_handler(E_ERR,'get_state_time','routine not written',source,revision,revdate)

    if ( .not. module_initialized ) call static_init_model
    time=cosmo_fc_time

    return

  end function get_state_time


!------------------------------------------------------------------------
!>


  function get_state_vector() result (sv)

    real(r8)             :: sv(1:model_size)

    integer              :: islab,ikind,nx,ny,sidx,eidx
    real(r8),allocatable :: mydata(:,:)

  call error_handler(E_ERR,'get_state_vector','routine not written',source,revision,revdate)

    if ( .not. module_initialized ) call static_init_model

    call set_allowed_state_vector_vars()

    do islab=1,nslabs
      ikind=cosmo_slabs(islab)%dart_kind
      if (ikind>0) then
        if (is_allowed_state_vector_var(ikind)) then
          nx=state_vector_vars(ikind)%nx
          ny=state_vector_vars(ikind)%ny
          allocate(mydata(1:nx,1:ny))
          mydata=get_data_from_binary(cosmo_restart_file,grib_header(islab),nx,ny)
          sidx=cosmo_slabs(islab)%dart_sindex
          eidx=cosmo_slabs(islab)%dart_eindex
          state_vector(sidx:eidx)=reshape(mydata,(/ (nx*ny) /))
          deallocate(mydata)
        endif
      endif
    enddo

    sv(:)=state_vector(:)

    return

  end function get_state_vector


!------------------------------------------------------------------------
!>


  subroutine write_grib_file(sv,nfile)

    real(r8),intent(in)           :: sv(:)
    character(len=128),intent(in) :: nfile

    integer                       :: istat = 0
    integer                       :: islab,ipos,lpos,bpos,griblen
    integer                       :: mylen,hlen,ix,iy,nx,ny,idx
    integer                       :: dval,ibsf,idsf
    real(r8)                      :: bsf,dsf
    integer(kind=1)               :: bin4(4)
    integer(kind=1),allocatable   :: bytearr(:),tmparr(:)
    real(r8),allocatable          :: mydata(:,:)
    real(r8)                      :: ref_value
    integer                       :: gribunitin, gribunitout,irec
    integer                       :: recpos(nslabs+1),bytpos(nslabs+1),myrlen,myblen

    logical                       :: write_from_sv = .false.

  call error_handler(E_ERR,'write_grib_file','routine not written',source,revision,revdate)
    if ( .not. module_initialized ) call static_init_model

    mylen=0

    ! read information on record and byte positions in GRIB file
    myrlen=(grib_header(2)%start_record-grib_header(1)%start_record)
    myblen=(grib_header(2)%start_record-grib_header(1)%start_record)*4+grib_header(2)%data_offset
    DO islab=1,nslabs
       recpos(islab)=grib_header(islab)%start_record
       bytpos(islab)=(grib_header(islab)%start_record-1)*4+1+grib_header(islab)%data_offset+4
    enddo
    recpos(nslabs+1)=recpos(nslabs)+myrlen
    bytpos(nslabs+1)=bytpos(nslabs)+myblen

    ! generate byte_array to write
    ALLOCATE(bytearr(1:bytpos(nslabs+1)+100))
    bytearr(:)=0

    ! read all data from the input GRIB file
    gribunitin  = get_unit()
    OPEN(gribunitin,FILE=TRIM(cosmo_restart_file),FORM='UNFORMATTED',ACCESS='DIRECT',RECL=record_length)
    irec=1
    lpos=1

    do while (istat==0)
      read(gribunitin,rec=irec,iostat=istat) bin4
      irec=irec+1
      lpos=lpos+4
      bytearr(lpos:lpos+3)=bin4
    enddo
    griblen=lpos-1

    call close_file(gribunitin)

    ipos=bytpos(1)

    if ( debug > 0 .and. do_output() ) write(*,*)'number of slabs is ',nslabs

    do islab=1,nslabs

      if ( debug > 0 .and. do_output() ) write(*,'(A8,A,I4,A,I4,A,I12)')cosmo_slabs(islab)%varname_short," is GRIB record ",islab," of ",nslabs,", byte position is ",ipos

      nx=cosmo_slabs(islab)%dims(1)
      ny=cosmo_slabs(islab)%dims(2)

      idx=cosmo_slabs(islab)%dart_sindex

      ! check if variable is in state vector

      write_from_sv = .false.
      if (idx >= 0) write_from_sv = .true.

      if (write_from_sv) then

        ! if variable is in state vector

        if ( debug > 0 .and. do_output() ) write(*,*)'         ... data is written to GRIB file from state vector'

        if ( debug > 0 .and. do_output() ) write(*,'(A8,A,I4,A,I4,A,2(1x,I12))')cosmo_slabs(islab)%varname_short," is GRIB record ",islab," of ",nslabs,", i1/i2 are ",idx,(idx+nx*ny-1)

        allocate(mydata(1:nx,1:ny))
        mydata(:,:)=reshape(sv(idx:(idx+nx*ny-1)),(/ nx,ny /))

        if (cosmo_slabs(islab)%dart_kind==KIND_PRESSURE_PERTURBATION) then
          mydata(:,:)=mydata(:,:)-non_state_data%p0fl(:,:,cosmo_slabs(islab)%ilevel)
        endif

        ref_value=minval(mydata)
        bin4(1:4)=from_float1(ref_value,cosmo_slabs(islab)%ref_value_char)

        hlen=size(grib_header(islab)%pds)+size(grib_header(islab)%gds)

        ! offset to binary data section in GRIB data, 8 because of indicator section
        bpos=ipos+hlen+8

        ! set new reference value in GRIB data
        bytearr(bpos+6:bpos+9)=bin4

        ! get the binary scale factor
        CALL byte_to_word_signed(bytearr(bpos+4:bpos+5),ibsf,2)
        bsf=FLOAT(ibsf)

        ! get the decimal scale factor
        CALL byte_to_word_signed(grib_header(islab)%pds(27:28),idsf,2)
        dsf=FLOAT(idsf)

        ! allocate a temporal array to save the binary data values
        allocate(tmparr((nx*ny*2)))
        lpos=1
        DO iy=1,ny
        DO ix=1,nx
          dval=int((mydata(ix,iy)-ref_value)*((10.**dsf)/(2.**bsf)))
          tmparr(lpos:lpos+1)=word_to_byte(dval)
          lpos=lpos+2
        enddo
        enddo

        deallocate(mydata)

        ! overwrite the old with new data
        bytearr(bpos+11:(bpos+nx*ny*2))=tmparr(1:(nx*ny*2))

        deallocate(tmparr)

        ipos=bytpos(islab+1)

      else

        if ( debug > 0 .and. do_output() ) write(*,*)'         ... data is copied from old grib file'

        ! if variable is not in state vector then skip this slab

        ipos=bytpos(islab+1)

      endif

    enddo

    ! write the new GRIB file
    gribunitout = get_unit()
    OPEN(gribunitout,FILE=TRIM(nfile),FORM='UNFORMATTED',ACCESS='stream')
    WRITE(gribunitout) bytearr(5:griblen)

    call close_file(gribunitout)

    return

  end subroutine write_grib_file


!------------------------------------------------------------------------
!>


function get_cosmo_filename(filetype)
character(len=*), optional, intent(in) :: filetype
character(len=256) :: get_cosmo_filename
character(len=256) :: lj_filename

call error_handler(E_ERR,'get_cosmo_filename','routine not written',source,revision,revdate)

lj_filename = adjustl(cosmo_restart_file)

if (present(filetype)) then
   if (trim(filetype) == 'netcdf') then
      lj_filename = adjustl(cosmo_netcdf_file)
   endif
endif

get_cosmo_filename = trim(lj_filename)

end function get_cosmo_filename


!------------------------------------------------------------------------
!>


  subroutine write_state_times(iunit, statetime, advancetime)
    integer,         intent(in) :: iunit
    type(time_type), intent(in) :: statetime, advancetime

    character(len=32) :: timestring
    integer           :: iyear, imonth, iday, ihour, imin, isec
    integer           :: ndays, nhours, nmins, nsecs
    type(time_type)   :: interval

  call error_handler(E_ERR,'write_state_times','routine not written',source,revision,revdate)

    call get_date(statetime, iyear, imonth, iday, ihour, imin, isec)
    write(timestring, "(I4,5(1X,I2))") iyear, imonth, iday, ihour, imin, isec
    write(iunit, "(A)") trim(timestring)

    call get_date(advancetime, iyear, imonth, iday, ihour, imin, isec)
    write(timestring, "(I4,5(1X,I2))") iyear, imonth, iday, ihour, imin, isec
    write(iunit, "(A)") trim(timestring)

    interval = advancetime - statetime
    call get_time(interval, nsecs, ndays)
    nhours = nsecs / (60*60)
    nsecs  = nsecs - (nhours * 60*60)
    nmins  = nsecs / 60
    nsecs  = nsecs - (nmins * 60)

    write(timestring, "(I4,3(1X,I2))") ndays, nhours, nmins, nsecs
    write(iunit, "(A)") trim(timestring)

  end subroutine write_state_times


!------------------------------------------------------------------------
!>


subroutine get_cosmo_grid(filename)

character(len=*), intent(in) :: filename

! Our little test case had these dimensions - only the names are important
!       time = UNLIMITED ; // (1 currently)
!       bnds = 2 ;
!       rlon = 30 ;
!       rlat = 20 ;
!       srlon = 30 ;
!       srlat = 20 ;
!       level = 50 ;
!       level1 = 51 ;
!       soil = 7 ;
!       soil1 = 8 ;

integer :: ncid, io, dimid

io = nf90_open(filename, NF90_NOWRITE, ncid)
call nc_check(io, 'get_cosmo_grid','open "'//trim(filename)//'"')

write(string1,*) 'time: '//trim(filename)
io = nf90_inq_dimid(ncid, 'time', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=ntime)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'bnds: '//trim(filename)
io = nf90_inq_dimid(ncid, 'bnds', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nbnds)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'rlon: '//trim(filename)
io = nf90_inq_dimid(ncid, 'rlon', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nrlon)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'rlat: '//trim(filename)
io = nf90_inq_dimid(ncid, 'rlat', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nrlat)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'srlon: '//trim(filename)
io = nf90_inq_dimid(ncid, 'srlon', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsrlon)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'srlat: '//trim(filename)
io = nf90_inq_dimid(ncid, 'srlat', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsrlat)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'level: '//trim(filename)
io = nf90_inq_dimid(ncid, 'level', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nlevel)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'level1: '//trim(filename)
io = nf90_inq_dimid(ncid, 'level1', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nlevel1)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'soil: '//trim(filename)
io = nf90_inq_dimid(ncid, 'soil', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsoil)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

write(string1,*) 'soil1: '//trim(filename)
io = nf90_inq_dimid(ncid, 'soil1', dimid)
call nc_check(io, 'get_cosmo_grid','inq_dimid '//trim(string1))
io = nf90_inquire_dimension(ncid, dimid, len=nsoil1)
call nc_check(io, 'get_cosmo_grid','inquire_dimension '//trim(string1))

if ( debug > 5 .and. do_output() ) then

   write(string1,*)'time   has dimension ',ntime
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'bnds   has dimension ',nbnds
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'rlon   has dimension ',nrlon
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'rlat   has dimension ',nrlat
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'srlon  has dimension ',nsrlon
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'srlat  has dimension ',nsrlat
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'level  has dimension ',nlevel
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'level1 has dimension ',nlevel1
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'soil   has dimension ',nsoil
   call error_handler(E_MSG,'get_cosmo_grid',string1)

   write(string1,*)'soil1  has dimension ',nsoil1
   call error_handler(E_MSG,'get_cosmo_grid',string1)

endif

call get_grid_var(ncid,  'lon' ,  nrlon,  nrlat, filename)
call get_grid_var(ncid,  'lat' ,  nrlon,  nrlat, filename)
call get_grid_var(ncid, 'slonu', nsrlon,  nrlat, filename)
call get_grid_var(ncid, 'slatu', nsrlon,  nrlat, filename)
call get_grid_var(ncid, 'slonv',  nrlon, nsrlat, filename)
call get_grid_var(ncid, 'slatv',  nrlon, nsrlat, filename)

where(lon <   0.0_r8) lon = lon + 360.0_r8
where(lat < -90.0_r8) lat = -90.0_r8
where(lat >  90.0_r8) lat =  90.0_r8

if ( debug > 5 .and. do_output() ) then
   ! call some sort of summary routine for min/max of grid vars
endif

call get_vcoord(ncid, filename)

call nc_check(nf90_close(ncid),'get_cosmo_grid','close "'//trim(filename)//'"' )

end subroutine get_cosmo_grid


!------------------------------------------------------------------------
!> all grid variables in the netCDF files are 2D
!> this does not handle scale, offset, missing_value, _FillValue etc.

subroutine get_grid_var(ncid,varstring,expected_dim1,expected_dim2,filename)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: varstring
integer,          intent(in) :: expected_dim1
integer,          intent(in) :: expected_dim2
character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
integer :: io, VarID, DimID, ndims
integer :: i

write(string3,*)trim(varstring)//' '//trim(filename)

call nc_check(nf90_inq_varid(ncid, trim(varstring), VarID), &
         'get_grid_var', 'inq_varid '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims), &
         'get_grid_var', 'inquire_variable '//trim(string3))

!> @TODO more informative error message
if (ndims /= 2) then
   call error_handler(E_ERR,'get_grid_var','wrong shape for '//string3, &
              source, revision, revdate)
endif

DimensionLoop : do i = 1,ndims

   write(string1,'(''inquire dimension'',i2,A)') i,trim(string3)

   call nc_check(nf90_inquire_dimension(ncid, dimIDs(i), len=dimlens(i)), &
                                       'get_grid_var', string1)

enddo DimensionLoop

! Check that the variable actual sizes match the expected sizes,

if (dimlens(1) .ne. expected_dim1 .or. &
    dimlens(2) .ne. expected_dim2 ) then
   write(string1,*)'expected dimension 1 to be ',expected_dim1, 'was', dimlens(1)
   write(string2,*)'expected dimension 2 to be ',expected_dim2, 'was', dimlens(2)
   call error_handler(E_ERR, 'get_grid_var', string1, &
              source, revision, revdate, text2=string2)
endif

select case (trim(varstring))
   case ("lon")
      allocate( lon(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, lon)
      call nc_check(io, 'get_grid_var', 'get_var '//trim(string3))

   case ("lat")
      allocate( lat(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, lat)
      call nc_check(io, 'get_grid_var', 'get_var '//trim(string3))

   case ("slonu")
      allocate( slonu(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slonu)
      call nc_check(io, 'get_grid_var', 'get_var '//trim(string3))

   case ("slatu")
      allocate( slatu(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slatu)
      call nc_check(io, 'get_grid_var', 'get_var '//trim(string3))

   case ("slonv")
      allocate( slonv(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slonv)
      call nc_check(io, 'get_grid_var', 'get_var '//trim(string3))

   case ("slatv")
      allocate( slatv(dimlens(1), dimlens(2)) )
      io = nf90_get_var(ncid, VarID, slatv)
      call nc_check(io, 'get_grid_var', 'get_var '//trim(string3))

   case default
      write(string1,*)'unsupported grid variable '
      call error_handler(E_ERR,'get_grid_var', string1, &
                 source, revision, revdate, text2=string3)

end select


end subroutine get_grid_var

!------------------------------------------------------------------------
!>

subroutine get_vcoord(ncid, filename)
integer,          intent(in) :: ncid
character(len=*), intent(in) :: filename

integer, dimension(NF90_MAX_VAR_DIMS) :: dimIDs
integer, dimension(NF90_MAX_VAR_DIMS) :: dimlens
integer :: io, VarID, DimID, ndims

write(string3,*)' vcoord from '//trim(filename)

call nc_check(nf90_inq_varid(ncid, 'vcoord', VarID), &
         'get_vcoord', 'inq_varid '//trim(string3))

call nc_check(nf90_inquire_variable(ncid, VarID, dimids=dimIDs, ndims=ndims), &
         'get_vcoord', 'inquire_variable '//trim(string3))

!> @TODO more informative error message
if (ndims /= 1) then
   call error_handler(E_ERR,'get_vcoord','wrong shape for '//string3, &
              source, revision, revdate)
endif

io = nf90_inquire_dimension(ncid, dimIDs(1), len=dimlens(1))
call nc_check(io, 'get_vcoord', 'inquire_dimension '//trim(string3))

! Check that the variable actual sizes match the expected sizes,

if (dimlens(1) .ne. nlevel1 ) then
   write(string1,*)'expected dimension to be ',nlevel1, 'was', dimlens(1)
   call error_handler(E_ERR, 'get_vcoord', string1, source, revision, revdate)
endif

allocate( vcoord%level1(nlevel1) )
io = nf90_get_var(ncid, VarID, vcoord%level1 )
call nc_check(io, 'get_vcoord', 'get_var '//trim(string3))

io = nf90_get_att(ncid, VarID, 'long_name', vcoord%long_name)
call nc_check(io, 'get_vcoord', 'get_att long_name '//trim(string3))

io = nf90_get_att(ncid, VarID, 'units', vcoord%units)
call nc_check(io, 'get_vcoord', 'get_att units '//trim(string3))

io = nf90_get_att(ncid, VarID, 'ivctype', vcoord%ivctype)
call nc_check(io, 'get_vcoord', 'get_att ivctype '//trim(string3))

io = nf90_get_att(ncid, VarID, 'irefatm', vcoord%irefatm)
call nc_check(io, 'get_vcoord', 'get_att irefatm '//trim(string3))

io = nf90_get_att(ncid, VarID, 'p0sl', vcoord%p0sl)
call nc_check(io, 'get_vcoord', 'get_att p0sl '//trim(string3))

io = nf90_get_att(ncid, VarID, 't0sl', vcoord%t0sl)
call nc_check(io, 'get_vcoord', 'get_att t0sl '//trim(string3))

io = nf90_get_att(ncid, VarID, 'dt0lp', vcoord%dt0lp)
call nc_check(io, 'get_vcoord', 'get_att dt0lp '//trim(string3))

io = nf90_get_att(ncid, VarID, 'vcflat', vcoord%vcflat)
call nc_check(io, 'get_vcoord', 'get_att vcflat '//trim(string3))

io = nf90_get_att(ncid, VarID, 'delta_t', vcoord%delta_t)
call nc_check(io, 'get_vcoord', 'get_att delta_t '//trim(string3))

io = nf90_get_att(ncid, VarID, 'h_scal', vcoord%h_scal)
call nc_check(io, 'get_vcoord', 'get_att h_scal '//trim(string3))

! Print a summary
if (debug > 5 .and. do_output()) then
   write(logfileunit,*)
   write(logfileunit,*)'vcoord long_name: ',trim(vcoord%long_name)
   write(logfileunit,*)'vcoord     units: ',vcoord%units
   write(logfileunit,*)'vcoord   ivctype: ',vcoord%ivctype
   write(logfileunit,*)'vcoord   irefatm: ',vcoord%irefatm
   write(logfileunit,*)'vcoord      p0sl: ',vcoord%p0sl
   write(logfileunit,*)'vcoord      t0sl: ',vcoord%t0sl
   write(logfileunit,*)'vcoord     dt0lp: ',vcoord%dt0lp
   write(logfileunit,*)'vcoord    vcflat: ',vcoord%vcflat
   write(logfileunit,*)'vcoord   delta_t: ',vcoord%delta_t
   write(logfileunit,*)'vcoord    h_scal: ',vcoord%h_scal

   write(*,*)
   write(*,*)'vcoord long_name: ',trim(vcoord%long_name)
   write(*,*)'vcoord     units: ',vcoord%units
   write(*,*)'vcoord   ivctype: ',vcoord%ivctype
   write(*,*)'vcoord   irefatm: ',vcoord%irefatm
   write(*,*)'vcoord      p0sl: ',vcoord%p0sl
   write(*,*)'vcoord      t0sl: ',vcoord%t0sl
   write(*,*)'vcoord     dt0lp: ',vcoord%dt0lp
   write(*,*)'vcoord    vcflat: ',vcoord%vcflat
   write(*,*)'vcoord   delta_t: ',vcoord%delta_t
   write(*,*)'vcoord    h_scal: ',vcoord%h_scal
endif

end subroutine get_vcoord


!------------------------------------------------------------------------

!>  This routine checks the user input against the variables available in the
!>  input netcdf file to see if it is possible to construct the DART state vector
!>  specified by the input.nml:model_nml:clm_variables  variable.
!>  Each variable must have 6 entries.
!>  1: GRIB table version number
!>  2: variable name
!>  3: DART KIND
!>  4: minimum value - as a character string - if none, use 'NA'
!>  5: maximum value - as a character string - if none, use 'NA'
!>  6: does the variable get updated in the restart file or not ...
!>     all variables will be updated INTERNALLY IN DART
!>     'UPDATE'       => update the variable in the restart file
!>     'NO_COPY_BACK' => do not copy the variable back to the restart file

subroutine parse_variable_table( state_variables, ngood, table )

character(len=*), dimension(:),   intent(in)  :: state_variables
integer,                          intent(out) :: ngood
character(len=*), dimension(:,:), intent(out) :: table

integer :: nrows, ncols, ivar
character(len=NF90_MAX_NAME) :: gribtableversion
character(len=NF90_MAX_NAME) :: gribvar
character(len=NF90_MAX_NAME) :: varname
character(len=NF90_MAX_NAME) :: dartstr
character(len=NF90_MAX_NAME) :: minvalstring
character(len=NF90_MAX_NAME) :: maxvalstring
character(len=NF90_MAX_NAME) :: state_or_aux

nrows = size(table,1)
ncols = size(table,2)

! This loop just repackages the 1D array of values into a 2D array.
! We can do some miniminal checking along the way.
! Determining which file to check is going to be more complicated.

ngood = 0
MyLoop : do ivar = 1, nrows

   gribtableversion = trim(state_variables(ncols*ivar - 6))
   gribvar          = trim(state_variables(ncols*ivar - 5))
   varname          = trim(state_variables(ncols*ivar - 4))
   dartstr          = trim(state_variables(ncols*ivar - 3))
   minvalstring     = trim(state_variables(ncols*ivar - 2))
   maxvalstring     = trim(state_variables(ncols*ivar - 1))
   state_or_aux     = trim(state_variables(ncols*ivar    ))

   call to_upper(state_or_aux)

   table(ivar,VT_GRIBVERSIONINDX) = trim(gribtableversion)
   table(ivar,VT_GRIBVARINDX)     = trim(gribvar)
   table(ivar,VT_VARNAMEINDX)     = trim(varname)
   table(ivar,VT_KINDINDX)        = trim(dartstr)
   table(ivar,VT_MINVALINDX)      = trim(minvalstring)
   table(ivar,VT_MAXVALINDX)      = trim(maxvalstring)
   table(ivar,VT_STATEINDX)       = trim(state_or_aux)

   ! If the first element is empty, we have found the end of the list.
   if ( table(ivar,1) == ' ' ) exit MyLoop

   ! Any other condition is an error.
   if ( any(table(ivar,:) == ' ') ) then
      string1 = 'input.nml &model_nml:variables not fully specified'
      string2 = 'must be 7 entries per variable. Last known variable name is'
      string3 = '['//trim(table(ivar,1))//'] ... (without the [], naturally)'
      call error_handler(E_ERR, 'parse_variable_table', string1, &
         source, revision, revdate, text2=string2, text3=string3)
   endif

   ! Make sure DART kind is valid

   if( get_raw_obs_kind_index(dartstr) < 0 ) then
      write(string1,'(''there is no obs_kind <'',a,''> in obs_kind_mod.f90'')') trim(dartstr)
      call error_handler(E_ERR,'parse_variable_table',string1,source,revision,revdate)
   endif

   ! Record the contents of the DART state vector

   if (debug > 8 .and. do_output()) then
      write(logfileunit,*)'variable ',ivar,' is ',trim(table(ivar,1)), ' ', trim(table(ivar,2)),' ', &
                                                  trim(table(ivar,3)), ' ', trim(table(ivar,4)),' ', &
                                                  trim(table(ivar,5)), ' ', trim(table(ivar,6)),' ', trim(table(ivar,7))
      write(     *     ,*)'variable ',ivar,' is ',trim(table(ivar,1)), ' ', trim(table(ivar,2)),' ', &
                                                  trim(table(ivar,3)), ' ', trim(table(ivar,4)),' ', &
                                                  trim(table(ivar,5)), ' ', trim(table(ivar,6)),' ', trim(table(ivar,7))
   endif

   ngood = ngood + 1

enddo MyLoop

if (ngood == nrows) then
   string1 = 'WARNING: There is a possibility you need to increase ''max_state_variables'''
   write(string2,'(''WARNING: you have specified at least '',i4,'' perhaps more.'')')ngood
   call error_handler(E_MSG,'parse_variable_table',string1,text2=string2)
endif

end subroutine parse_variable_table


!------------------------------------------------------------------------
!>

!> SetVariableAttributes() converts the information in the variable_table
!> to the progvar structure for each variable.
!> If the numerical limit does not apply, it is set to MISSING_R8, even if
!> it is the maximum that does not apply.

subroutine SetVariableAttributes(ivar)

integer, intent(in) :: ivar

integer  :: ios, ivalue
real(r8) :: minvalue, maxvalue

progvar(ivar)%varname     = trim(variable_table(ivar,VT_VARNAMEINDX))
progvar(ivar)%kind_string = trim(variable_table(ivar,VT_KINDINDX))
progvar(ivar)%dart_kind   = get_raw_obs_kind_index( progvar(ivar)%kind_string )
progvar(ivar)%maxlevels   = 0
progvar(ivar)%dimlens     = 0
progvar(ivar)%dimnames    = ' '
progvar(ivar)%rangeRestricted   = BOUNDED_NONE
progvar(ivar)%minvalue          = MISSING_R8
progvar(ivar)%maxvalue          = MISSING_R8
progvar(ivar)%update            = .false.

if (variable_table(ivar,VT_STATEINDX)  == 'UPDATE') progvar(ivar)%update = .true.

! set the default values

minvalue = MISSING_R8
maxvalue = MISSING_R8
progvar(ivar)%minvalue = MISSING_R8
progvar(ivar)%maxvalue = MISSING_R8

! If the character string can be interpreted as an r8, great.
! If not, there is no value to be used.

read(variable_table(ivar,VT_MINVALINDX),*,iostat=ios) minvalue
if (ios == 0) progvar(ivar)%minvalue = minvalue

read(variable_table(ivar,VT_MAXVALINDX),*,iostat=ios) maxvalue
if (ios == 0) progvar(ivar)%maxvalue = maxvalue

read(variable_table(ivar,VT_GRIBVERSIONINDX),*,iostat=ios) ivalue
if (ios == 0) progvar(ivar)%tableID = ivalue

read(variable_table(ivar,VT_GRIBVARINDX),*,iostat=ios) ivalue
if (ios == 0) progvar(ivar)%variableID = ivalue

! rangeRestricted == BOUNDED_NONE  == 0 ... unlimited range
! rangeRestricted == BOUNDED_BELOW == 1 ... minimum, but no maximum
! rangeRestricted == BOUNDED_ABOVE == 2 ... maximum, but no minimum
! rangeRestricted == BOUNDED_BOTH  == 3 ... minimum and maximum

if (   (progvar(ivar)%minvalue /= MISSING_R8) .and. &
       (progvar(ivar)%maxvalue /= MISSING_R8) ) then
   progvar(ivar)%rangeRestricted = BOUNDED_BOTH

elseif (progvar(ivar)%maxvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_ABOVE

elseif (progvar(ivar)%minvalue /= MISSING_R8) then
   progvar(ivar)%rangeRestricted = BOUNDED_BELOW

else
   progvar(ivar)%rangeRestricted = BOUNDED_NONE

endif

! Check to make sure min is less than max if both are specified.

if ( progvar(ivar)%rangeRestricted == BOUNDED_BOTH ) then
   if (maxvalue < minvalue) then
      write(string1,*)'&model_nml state_variable input error for ',trim(progvar(ivar)%varname)
      write(string2,*)'minimum value (',minvalue,') must be less than '
      write(string3,*)'maximum value (',maxvalue,')'
      call error_handler(E_ERR,'SetVariableAttributes',string1, &
         source,revision,revdate,text2=trim(string2),text3=trim(string3))
   endif
endif

end subroutine SetVariableAttributes



!------------------------------------------------------------------------
!>


subroutine set_variable_binary_properties(dartid)
integer, intent(in) :: dartid

! Table to decode the record contents of igdsbuf
integer, parameter :: indx_numEW    =  5, &
                      indx_numNS    =  6, &
                      indx_startlat =  7, &
                      indx_startlon =  8, &
                      indx_endlat   = 10, &
                      indx_endlon   = 11

! Table to decode the record contents of ipdsbuf
!> @TODO no seconds?
integer, parameter :: indx_gribver   =   2, &
                      indx_var       =   7, &
                      indx_zlevtyp   =   8, &
                      indx_zlevtop   =   9, &
                      indx_zlevbot   =  10, &
                      indx_year      =  11, &
                      indx_mm        =  12, &
                      indx_dd        =  13, &
                      indx_hh        =  14, &
                      indx_min       =  15, &
                      indx_startstep =  17, &
                      indx_endstep   =  18, &
                      indx_nztri     =  19, &
                      indx_cc        =  22

integer, parameter :: NPDS  = 321   ! dimension of product definition section
integer, parameter :: NGDS  = 626   ! dimension of grid description section

integer(i4) :: ipdsbuf(NPDS)        ! pds: product definition section
integer(i4) :: igdsbuf(NGDS)        ! gds: grid definition section

integer(i4) :: nudat,      &        ! unit number
               izerr,      &        ! error status
               ivar,       &        ! variable reference number based on iver
               iver,       &        ! version number of GRIB1 indicator table
               iz_countl            ! counter for binary data

real(r8) :: rbuf(nrlon*nrlat)       ! data to be read
real(r8) :: zvc_params(nlevel1)     ! height levels

real(r8) :: psm0,             &     ! initial value for mean surface pressure ps
            dsem0,            &     ! initial value for mean dry static energy
            msem0,            &     ! initial value for mean moist static energy
            kem0,             &     ! initial value for mean kinetic energy
            qcm0,             &     ! initial value for mean cloudwater content
            refatm_p0sl,      &     ! constant reference pressure on sea-level
            refatm_t0sl,      &     ! constant reference temperature on sea-level
            refatm_dt0lp,     &     ! d (t0) / d (ln p0)
            refatm_delta_t,   &     ! temperature difference between sea level and stratosphere (for irefatm=2)
            refatm_h_scal,    &     ! scale height (for irefatm=2)
            refatm_bvref,     &     ! constant Brund-Vaisala-frequency for irefatm=3
            vcoord_vcflat           ! coordinate where levels become flat

integer(i4) :: ntke,          &     ! time level for TKE
               izvctype_read        ! check vertical coordinate type in restarts

integer  :: fid, i, nx, ny, ilev, ilevp1, ilevtyp
integer  :: icc, iyy, imm, idd, ihh, imin, iccyy
integer  :: istartstep, iendstep, nztri
real(r8) :: lat1, latN, lon1, lonN

! variables that record what kind of variable is being read
! must switch on change of variable 
integer :: old_tableID, old_varID, old_izvctype

!-----------------------------------------------------------------------

fid = open_file(cosmo_restart_file, form='unformatted', action='read')

read(fid, iostat=izerr) psm0, dsem0, msem0, kem0, qcm0, ntke
if (izerr /= 0) then
   write(string1,*)'unable to read first record while searching for '//trim(progvar(dartid)%varname)
   call error_handler(E_ERR,'set_variable_binary_properties', string1, source, revision, revdate)
endif

read(fid, iostat=izerr) izvctype_read, refatm_p0sl, refatm_t0sl,   &
                        refatm_dt0lp, vcoord_vcflat, zvc_params
if (izerr /= 0) then
   write(string1,*)'unable to read second record while searching for '//trim(progvar(dartid)%varname)
   call error_handler(E_ERR,'set_variable_binary_properties', string1, source, revision, revdate)
endif

if     ( izvctype_read > 0 .and. izvctype_read <= 100 ) then
   !write(*,*) "izvctype_read = ", izvctype_read
   continue

elseif ( (izvctype_read > 100) .and. (izvctype_read <= 200) ) then

   read(fid, iostat=izerr) refatm_delta_t, refatm_h_scal
   if (izerr /= 0) then
      write(string1,*)'izvctype_read is ',izvctype_read, 'requiring us to read "refatm_delta_t, refatm_h_scal"'
      write(string2,*)'unable to read record while searching for '//trim(progvar(dartid)%varname)
      call error_handler(E_ERR,'set_variable_binary_properties', string1, source, revision, revdate, text2=string2)
   endif

elseif ( (izvctype_read > 200) .and. (izvctype_read <= 300) ) then

   read(fid, iostat=izerr) refatm_bvref
   if (izerr /= 0) then
      write(string1,*)'izvctype_read is ',izvctype_read, 'requiring us to read "refatm_bvref"'
      write(string2,*)'unable to read record while searching for '//trim(progvar(dartid)%varname)
      call error_handler(E_ERR,'set_variable_binary_properties', string1, source, revision, revdate, text2=string2)
   endif

else
   write(string1,*) 'izvctype_read is ',izvctype_read,' is unsupported.'
   call error_handler(E_ERR,'set_variable_binary_properties', string1, source, revision, revdate)
endif

!------------------------------------------------------------------------------
!Section 3: READ ALL RECORDS
!------------------------------------------------------------------------------

iz_countl    = 0
old_tableID  = 0
old_varID    = 0
old_izvctype = 0

read_loop: DO

   read(fid, iostat=izerr) ipdsbuf, igdsbuf, rbuf

   if (izerr < 0) then
     exit read_loop
   elseif (izerr > 0) then
      write(string1,*) 'ERROR READING RESTART FILE around data record ',iz_countl
      call error_handler(E_ERR,'set_variable_binary_properties', string1, source, revision, revdate)
   endif

   iz_countl = iz_countl + 1

   ! Decode ipdsbuf

   iver    = ipdsbuf(indx_gribver)
   ivar    = ipdsbuf(indx_var)
   ilevtyp = ipdsbuf(indx_zlevtyp)
   ilev    = ipdsbuf(indx_zlevtop)
   ilevp1  = ipdsbuf(indx_zlevbot)

   icc     = ipdsbuf(indx_cc)-1
   iyy     = ipdsbuf(indx_year)
   iccyy   = iyy + icc*100
   imm     = ipdsbuf(indx_mm)
   idd     = ipdsbuf(indx_dd)
   ihh     = ipdsbuf(indx_hh)
   imin    = ipdsbuf(indx_min)

   istartstep = ipdsbuf(indx_startstep)
   iendstep   = ipdsbuf(indx_endstep)
   nztri      = ipdsbuf(indx_nztri)

   write(*,'(A,3x,i8,5(1x,i2),4(1x,i4))') 'time for iz_countl = ', &
         iz_countl,icc,iyy,imm,idd,ihh,imin,nztri,istartstep,iendstep

   ! Decode igdsbuf

   nx   = igdsbuf(indx_numEW)
   ny   = igdsbuf(indx_numNS)
   lat1 = real(igdsbuf(indx_startlat),r8) * 0.001_r8
   lon1 = real(igdsbuf(indx_startlon),r8) * 0.001_r8
   latN = real(igdsbuf(indx_endlat  ),r8) * 0.001_r8
   lonN = real(igdsbuf(indx_endlon  ),r8) * 0.001_r8

   ! Since Fortran binary reads have no way to know if they fail or not
   ! we are going to compare our input with the sizes encoded in the file.

   if (nx .ne. nrlon .or. ny .ne. nrlat) then
      write(string1,*) 'reading variable ',trim(progvar(dartid)%varname), ' from '//trim(cosmo_restart_file)
      write(string2,*) '(file) nx /= nrlon ',nx,nrlon, ' or '
      write(string3,*) '(file) ny /= nrlat ',ny,nrlat
      call error_handler(E_ERR,'set_variable_binary_properties', string1, &
                 source, revision, revdate, text2=string2, text3=string3)
   endif

   ! see if the slab we are reading defines the start of a new variable

   if (old_tableID  .ne. iver .or. &
       old_varID    .ne. ivar .or. &
       old_izvctype .ne. ilevtyp) then
      write (*,'(A,2x,8(1x,i4),2(1xf15.8))')'slab ',iz_countl, iver, ivar, nx, ny, &
                ilev, ilevp1, ilevtyp, MINVAL(rbuf), MAXVAL(rbuf)
      ! Check to see if the slab/variable/etc. is one we want
      ! if so - count up the slabs, add to the variable size and keep moving
   endif

   old_tableID  = iver
   old_varID    = ivar
   old_izvctype = ilevtyp

enddo read_loop

close(fid, IOSTAT=izerr)
if (izerr /= 0) then
   write(string1,*) 'closing '//trim(cosmo_restart_file)//' while looking for ',trim(progvar(dartid)%varname)
   call error_handler(E_ERR,'set_variable_binary_properties', string1, source, revision, revdate)
endif

end subroutine set_variable_binary_properties


!------------------------------------------------------------------------
!>

subroutine progvar_summary

integer :: ivar,i

! only have 1 task write the report
if ( .not. do_output() ) return

do ivar = 1,nfields
   write(*,*)
   write(*,*)'variable ',ivar,' is ',trim(progvar(ivar)%varname)
   write(*,*)'   long_name       ',  trim(progvar(ivar)%long_name)
   write(*,*)'   units           ',  trim(progvar(ivar)%units)
   write(*,*)'   tableID         ',       progvar(ivar)%tableID
   write(*,*)'   variableID      ',       progvar(ivar)%variableID
   write(*,*)'   levtypeID       ',       progvar(ivar)%levtypeID
   write(*,*)'   maxlevels       ',       progvar(ivar)%maxlevels
   write(*,*)'   varsize         ',       progvar(ivar)%varsize
   write(*,*)'   index1          ',       progvar(ivar)%index1
   write(*,*)'   indexN          ',       progvar(ivar)%indexN
   write(*,*)'   dart_kind       ',       progvar(ivar)%dart_kind
   write(*,*)'   rangeRestricted ',       progvar(ivar)%rangeRestricted
   write(*,*)'   minvalue        ',       progvar(ivar)%minvalue
   write(*,*)'   maxvalue        ',       progvar(ivar)%maxvalue
   write(*,*)'   kind_string     ',trim(  progvar(ivar)%kind_string)
   write(*,*)'   update          ',       progvar(ivar)%update
   write(*,*)'   numdims         ',       progvar(ivar)%numdims
   do i = 1,progvar(ivar)%numdims
      write(*,*)'   dimension ',i,' has length ',progvar(ivar)%dimlens(i),' ',trim(progvar(ivar)%dimnames(i))
   enddo
enddo

end subroutine progvar_summary

!------------------------------------------------------------------------
!>

end module model_mod

! <next few lines under version control, do not edit>
! $URL: $
! $Id: $
! $Revision: $
! $Date: $
