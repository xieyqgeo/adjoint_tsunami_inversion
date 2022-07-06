
module VariableDefination



type :: GlobalParameters

    !************ Read From Config ************!
    real*8      ::  Version = 3.0 ! pcomcot version
    integer*4   ::  CoordinatesType        ! 0/spherical; 1/cartisian
    integer*4   ::  GoverningEquationType  ! 0/linear; 1/nonlinear; **only use linear
    real*8      ::  WaterDepthLimit        ! considered to be 0 if water depth is less than this
    integer*4   ::  InitialConditionType   ! Initial condition  0/file; 1/fault
    integer*4   ::  BoundaryConditionType  ! Boundary condition 0/open; 1/wall; **only use wall boundaries
    real*8      ::  TotalTime              ! Total run time
    real*8      ::  dt                     ! Time step
    real*8      ::  DTSaveData             ! Time interval to save snapshots
    integer*4   ::  SaveFlux               ! 0/no; 1/save flux;
    integer*4   ::  ComputeGreen           ! 0/Normal Simulation; 1/Computing Green Functions
    integer*4   ::  MinGridsPerNode        ! Min grids on each computing node
    integer*4   ::  FeedbackToParentLayer  ! feedback or not
    
    !************ Computed for Convinience *************!
    integer*4         ::  nCalculations              ! Value set in "readConfig", "readFaultParameters"
    integer*4         ::  TotalTimeSteps             ! Value set in "readConfig"
    integer*4         ::  NDTSaveData                ! Value set in "readConfig"
    integer*4         ::  NumLayers                  ! Value set in "checkFiles"
    integer*4         ::  NumLayerLevels             ! Value set in "determineLayerDependency"
    integer*4         ::  TopLayer                   ! Value set in "determineLayerDependency"
    character(999)    ::  InitialConditionFileName   ! Value set in "checkFiles"
    integer*4         ::  InitialConditionFileFormat ! Value set in "checkFiles" 1/nf; 2/xyz
    integer*4         ::  StartEastWest              ! Vaule set in "getBathymetryDataSize";
                                                     ! 0/Cartisian; 1/spherical long:160E~-160W(200E); 2/special long: 340W(-20E)~20E
    integer*4         ::  NumStations                ! Value set in "readStations"
    integer*4         ::  NumFaults                  ! Value set in "readFaultParameters"
    real*8,dimension(:),pointer  ::  t               ! Value set when solving equations

    !************ Default System Parameters ************!
    character(999)    ::  COMCOTParametersFileName  = "pcomcot.ctl"
    character(999)    ::  BathymetryFilePrefix = "layer"                  ! InitialCondition: .nf/.xyz
    character(999)    ::  InitialConditionFilePrefix = "InitialElevation" ! InitialCondition: .nf/.xyz
    character(999)    ::  FaultParametersFileName = "FaultParameters.ctl" ! Fault parameters
    character(999)    ::  StationFileName = "Stations.ctl"                ! File for stations coordinates
    integer*4         ::  ComputeDivisionOpt = 1                          ! 1: divide domain averagely; 2: child layer on 1 node
    integer*4         ::  nRowBoundary = 1                                ! Sync boundary rows between subdomains
    integer*4         ::  nRowBoundaryFlux = 1                            ! Sync boundary rows between subdomains for flux
    integer*4         ::  nRowBathymetry = 3                              ! At least 3 cells > WaterDepthLimit
    integer*4         ::  nRowBathymetryBoundary = 10                     ! At least 10 cells > WaterDepthLimit from boundary
    real*8            ::  GRAV = 9.807
    real*8            ::  PI = 3.141593
    real*8            ::  R_EARTH = 6378000.0
    real*8            ::  OMEGA = 7.2921159E-5

    !************  MPI Variables  ************!
    integer*4         ::  irank, nsizeTotal, master
    integer*4         ::  MaxNX, MinNX, MaxNY, MinNY   ! Value set in "partitionDomain"; Max/Min grids on one node
    real*8            ::  CPUTimeInitial               ! Value set in "pcomcot"
    real*8            ::  CPUTimeInitialCalculation    ! Value set in "pcomcot"
    real*8            ::  CPUTime(999,5)               ! Value set in all subroutines; iNode;AllTime/Prepare/CalculationTime/Communication/WriteToDisk


end type GlobalParameters



type  ::  LayerParameters

    integer*4       ::  Level                  ! Value set in "determineLayerDependency"
    integer*4       ::  Parent                 ! Value set in "determineLayerDependency"
    character(999)  ::  BathymetryFileName     ! Value set in "getBathymetryDataSize"
    integer*4       ::  BathymetryFileFormat   ! Value set in "getBathymetryDataSize"
    real*8          ::  xmin, dx, xmax         ! Value set in "getBathymetryDataSize"
    real*8          ::  ymin, dy, ymax         ! Value set in "getBathymetryDataSize"
    integer*4       ::  NX, NY                 ! Value set in "getBathymetryDataSize"
    real*8, dimension(:), pointer     ::  X, Y ! Value set in "readBathymetry"
    real*8, dimension(:,:), pointer   ::  Z    ! Value set in "readBathymetry"
    real*8          ::  zmin, zmax             ! Value set in "readBathymetry"
    real*8          ::  dt                     ! Value set in "cflCheck"
    integer*4       ::  nStepsPerTimeStep      ! Value set in "cflCheck"
    real*8          ::  dtratio                ! Value set in "cflCheck"
    integer*4       ::  nsize                          ! Value set in "partitionDomain"; ComputeNodes for layer
    integer*4       ::  npartx, nparty                 ! Value set in "partitionDomain"; partition
    integer*4       ::  MaxNX, MinNX, MaxNY, MinNY     ! Value set in "partitionDomain"; Max/Min grids on one node
    integer*4       ::  PartitionInfo(999,4)           ! Value set in "partitionDomain"; ComputeNode; nstartx,nendx,nstarty,nendym
    integer*4       ::  BoundarySendRecvCount          ! Value set in "partitionDomain"
    integer*4       ::  BoundarySendRecv(999,14)       ! Value set in "partitionDomain"; Count; FromNode,ToNode,(nstartx,nendx...)*4 for H,M,N
    integer*4       ::  ParentToChildSendRecvCount     ! Value set in "partitionDomain";
    integer*4       ::  ParentToChildSendRecv(9999,5)  ! Value set in "partitionDomain"; Count; FromNode,ToNode,Flag(hxmin,hxmax,hymin,hymax),start,end
    integer*4       ::  ChildToParentSendRecvCount     ! Value set in "partitionDomain"
    integer*4       ::  ChildToParentSendRecv(9999,7)  ! Value set in "partitionDomain"; Count; FromNode,ToNode,H/M/N,(nstartx,nendx...)

    !**** Values Computed as Parameters on each node  ****!
    !**** Values set in "computeParameters" ****!
    real*8, dimension(:), pointer      ::  CPX, CPY  ! Coriolis parameters CP=2*Omega*SIN(Y)
    real*8                             ::  CC1, CC2, CC3, CC4     ! coefficients for Cartisian coordinates
    real*8, dimension(:), pointer      ::  CSY, CS1, CS2, CS3     ! coefficients for spherical coordinates
    real*8                             ::  CS4

    !**** Values Computed by Solving Equations  ****!
    real*8, dimension(:,:,:), pointer  ::  H   ! water elevation, value in last/this step; x; y position
    real*8, dimension(:,:,:), pointer  ::  M   ! volume flux in x direction (h*u)
    real*8, dimension(:,:,:), pointer  ::  N   ! volume flux in y direction (h*v)
    real*8, dimension(:,:), pointer    ::  hxmin, hxmax, hymin, hymax ! Temperoral boundary values; Value set in "getLayerBoundaryValueFromParent"

end type LayerParameters



type  ::  StationParameters

    real*8      ::  X,Y      ! Value set in "readStations";    Stations coordinates
    integer*4   ::  nLayer   ! Value set in "readStations";    In which layer is the station located
    integer*4   ::  nNode    ! Value set in "partitionDomain"; On which node is the station located

    real*8, dimension(:), pointer      ::  H, M, N  !Values at this station
    
end type StationParameters



type FaultParameters

    real*8         ::  T0           ! rupture starting time
    integer*4      ::  NT           ! time step when rupture starts
    real*8         ::  Depth        ! focal depth
    real*8         ::  Length       ! length of fault plane
    real*8         ::  Width        ! width of fault plane
    real*8         ::  Slip         ! dislocation
    real*8         ::  Rake         ! rake
    real*8         ::  HSlip        ! horizontal dislocation     HSlip = Slip*cos(Rake)
    real*8         ::  PSlip        ! perpendicular dislocation  PSlip = Slip*sin(Rake)
    real*8         ::  Strike       ! strike
    real*8         ::  Dip          ! dip
    real*8         ::  Y0           ! epicenter(latitude)
    real*8         ::  X0           ! epicenter(longitude)

end type FaultParameters



end module VariableDefination
