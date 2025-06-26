Coupled Application
===================

This is the code for a simple application in which two "models" are coupled
using the OASIS3-MCT coupler. A single field is passed from one model to the
other. In the example provded with this application, a field is passed from a
global toy LFRic model to a limited-area toy LFRic model. The regridding weights
for that particular coupling are included in the example directory. The models
run for a single timestep. There is no science in the model.

The same executable is used for both models. An optional command line argument
(```-c```) is used to configure the component name - which is then used to
determine whether the executable is sending or receiving the field data.

To run the implementation from the `example` directory, use a command like such
as:
```
mpiexec -np 1 ../bin/coupled configuration_glo.nml -c lfric_o : -np 1 ../bin/coupled configuration_lam.nml -c lfric_i
```

The application is based on the skeleton application with the following changes:

1. ```applications/coupled/source/algorithm/coupled_alg_mod.X90``` - this file
has been extensively changed. The original code (that created a field of
constant value and calculated the divergence of it - and checked it was zero)
has been replaced by iterations over the coupling send and receive field
collections, calling the coupler though the external field object for all
fields. (Note: this has changed from a ```.x90``` to a ```.X90```, as there are
now ```ifdef```s on the ```MCT``` preprocessor directive).
2. ```applications/coupled/source/driver/init_coupled_mod.X90``` - this now
initialises two fields: ```field_1``` for sending to the coupler and
```field_2``` for receiving from the coupler. The fields are added to the
depository and also the appropriate coupling field collection. ```field_1```
is initialised with the longitude (in degrees) of its points. ```field_2``` is
initialised with the latitude - but this is never used. (Note: this has changed
from a ```.F90``` to a ```.X90```, as the latitude/longitude of the data points
are generated from a kernel invoke).
3. ```applications/coupled/source/driver/coupled_driver_mod.f90``` - this is
updated so all fields will be 2d (to make coupling easier) and tidied the
interface to the algorithm by passing the whole of ```modeldb```, rather than
individual fields.
4. ```applications/coupled/source/coupled.f90``` - this has extra code to deal
with the new ```-c``` optional command line argument.
5. ```applications/coupled/example``` now contains:
    * ```configuration_glo.nml``` - holds the configuration namelists for the
global coupled skeleton model
    * ```configuration_lam.nml``` - holds the configuration namelists for the
limited-area coupled skeleton model
    * ```namcouple``` - holds the OASIS3-MCT configuration that describes the
coupling
    * ```mesh_C12.nc``` - NetCDF file that holds a single global C12
cubed-sphere mesh in UGRID format
    * ```mesh_LAM50x50-2x2.nc``` - NetCDF file that holds a single 50x50
limited-area planar mesh in UGRID format
    * ```rmp_lfric_glo_to_lfric_lam.nc```- NetCDF file that holds the remapping
weight data for coupling fields held on the global mesh to fields on the
limited-area mesh

```
Call tree for Coupled Application
=================================
(Just the calls involved in coupling - the rest is the same as the skeleton app)


coupled                                      (applications/coupled/source/coupled.f90)
|   (^main program unit)
|__init_comm                                 (components/driver/source/driver_comm_mod.F90)
|  |__cpl_ptr%initialise                     (components/coupling/source/coupling_mod.F90)
|          (^also sets up the MPI communicator)
|
|__initialise                                (applications/coupled/source/driver/coupled_driver_mod.f90)
|  |__init_coupled                           (applications/coupled/source/driver/init_coupled_mod.X90)
|     |    (^creates and initialises fields)
|     |__get_coupling_fields                 (components/coupling/source/coupling_mod.F90)
|     |       (^gets the list of fields from the namcouple)
|     |__(set up the cpl_snd and cpl_rcv field_collections)
|     |__cpl_ptr%define_partitions           (components/coupling/source/coupling_mod.F90)
|     |__cpl_ptr%define_variables            (components/coupling/source/coupling_mod.F90)
|     |__cpl_ptr%end_definition              (components/coupling/source/coupling_mod.F90)
|
|__step                                      (applications/coupled/source/driver/coupled_driver_mod.f90)
|  |__coupled_alg                            (applications/coupled/source/algorithm/coupled_alg_mod.F90)
|     |__(top of iteration loop over cpl_snd_2d or cpl_rcv_2d)
|     |__coupler_exchange_2d%initialise      (components/coupling/source/coupler_exchange_2d_mod.F90)
|     |__coupler_exchange_2d%set_time        (components/coupling/source/coupler_exchange_2d_mod.F90)
|     |__coupler_exchange_2d%copy_from_lfric (components/coupling/source/coupler_exchange_2d_mod.F90)
|     |                  (or copy_to_lfric)  (components/coupling/source/coupler_exchange_2d_mod.F90)
|     |__coupler_exchange_2d%clear           (components/coupling/source/coupler_exchange_2d_mod.F90)
|     |__(loop back to top of iteration loop)
|
|__finalise                                  (applications/coupled/source/driver/coupled_driver_mod.f90)
```
