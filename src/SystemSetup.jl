module SystemSetup

#using #dependencies

#import #methods to overload

export Bus, ConductorModel, ThermalModel, Line, Gen,
        Load, Wind, Storage, Params, Control, 
        OptParams, CostFunction, Optimize, Time, System, 
        mat2sys

# BY JONAS KERSULIS, Jan. 22 2015
# This module defines the types used to characterize a transmission system.

# FROM MADS'S WORK
# See runLevel3vFull.m.
# The majority of System is created in setup files which then call
# runLevel3vFull.m (e.g. testLevel3 which calls systemSetup_RTS96 and
# then calls runLevel3vFull).
# Fields which are marked with an asterisk '*' are initialized or added to
# System in runLevel3vFull.  Unmarked fields should be initialized prior to
# runLevel3vFull being called.

type Bus # (number of buses in network)x1
    num::Array{Integer} # bus number
    ty::Array{Integer}  # bus type (1 = PQ, 2 = PV, 3 = slack)
    vBase               # base voltage in kV
    vMag                # voltage magnitude in puu
    vAng                # voltage angle in rad
    g                   # shunt conductance in pu
    b                   # shunt susceptance in pu
end

type ConductorModel # 1x1 structure containing...
    R_unit              # unit resistance of conductor [ohms/m]
    I_max               # current limit of network connection [A] (note that this is related to sLim field above. sLim = 0 or very large for voltage base will cause I_max to have a VERY large current rating (~50K) in order to approximate infinity.
    D                   # diameter of complete conductor [mm]
    Al_m                # unit mass of aluminum portion of ACSR [kg/m]
    St_m                # unit mass of steel portion of ACSR [kg/m]
    bundle              # number of conductors used in the network connection

    ConductorModel() = new()
end

type ThermalModel # 1x1 structure containing...
    T_amb               # ambient temperature [C]
    T_max               # ampacity-limited steady-state temperature [C] (related to sLim and I_max)
    η_c                 # conductive heat factor (varies w/ temp) [W/m-C]
    η_r                 # radiative heat factor [W/m-C^4]
    η_s                 # linear solar heat gain rate [W/m]
    mCp                 # line heat capacity per unit length (J/m-C)
    τ                   # (unitless) coefficient for temperature dynamics
    ρ                   # (C-m/W) coefficient for temperature dynamics
    γ                   # (unitless) coefficient for temperature dynamics
    fLossLim            # per-phase active power loss linearization point [W/m]

    ThermalModel() = new()
end

type Line # (number of lines in network)x1 structure
    # Unified branch model is chosen as:
    #         k  tmk:1    p     Zkm      q               m 
    #         |----8 8------O-----ZZZ------O-------------|
    #   Tap   +             |              |             +   Fixed
    #   bus   Vk     ykm^sh Z              Z ykm^sh      Vm   bus
    #         -             |              |             -
    #        ---------------------------------------------
    # with tmk = tap*exp(j*phase) is pu (complex) tap-ratio
    # o when tmk = 1          ===> transmission line
    # o when tmk = real       ===> TCT ===> tap-ratio (~1) [assumes ykm^sh=0]
    # o when tmk = complex    ===> PST ===> angle (rads)   [assumes ykm^sh=0]

    from::Array{Integer}            # from bus number
    to::Array{Integer}              # to bus number
    R                               # resistance (pu)
    X                               # reactance (pu)
    B                               # shunt susceptance (pu)
    sLim                            # 3-phase power rating of line (MVA)
    length                          # length of line (miles) (Level 3 can estimate lengths if data not available)
    f                               # active power flow on line (pu)
    fLoss                           # active power losses on line (pu)
    q                               # reactive power flow on line (pu)
    qLoss                           # reactive power losses on line (pu)
    delT                            # degrees above thermal line limit (C)
    isXfrmr::Array{Bool}            # logical specifying if line is a transformer 
    isTCXfrmr::Array{Bool}          # logical specifying if transformer is tap-changing
    isPSXfrmr::Array{Bool}          # logical specifying if transformer is phase-shifting
    tap                             # tap ratio (for Xfrmrs)
    phase                           # phase ratio (for Xfrmrs) (rad)
    maxPhase                        # maximum phase shift (rad) 
    conductorModel::ConductorModel  # 1x1 structure containing...
    thermalModel::ThermalModel      # 1x1 structure containing...
    Δθ                              # width of PWL segments approximating abs(theta_{ij})
end

type Gen # (number of generators in network)x1 structure
    bus::Array{Integer}     # bus number of generator
    ID::Array{Integer}      # generator identifier (not needed)
    status::Array{Integer}  # 0 = out of service, 1 = in service (not needed: Level 3 assumes all gens are available)
    Pinj                    # active power injection (pu)
    Qinj                    # reactive power injection (pu)
    Pramp                   # active power ramping rate (pu/hr)
    vSet                    # voltage setpoint of generator (pu)
    Pmin                    # min active power injection (pu)
    Pmax                    # max active power injection (pu)
    Qmin                    # min reactive power injection (pu) (should be <= 0)
    Qmax                    # max reactive power injection (pu) (should be >= 0)
    pUpLim                  # active power ramp rate limit (pu/hr) (for increasing output)
    pDownLim                # active power ramp rate limit (pu/hr) (for decreasing output)
    qUpLim                  # reactive power ramp rate limit (pu/hr) (for increasing output)
    qDownLim                # reactive power ramp rate limit (pu/hr) (for decreasing output)
end

 type Load # (number of loads in network)x1 structure
    bus::Array{Integer}     # bus number of load
    Pload                   # active power demand forecast (pu) 1x(number of time intervals) vector
    Qload                   # reactive power demand (pu) 1x(number of time intervals) vector
    Pdemand                 # present active power load (pu)
    Qdemand                 # present reactive power load (pu)
    Pshed                   # active power shed from load through demand response (pu)
    Qshed                   # reactive power shed from load through demand response (pu)
end

type Wind # (number of wind farms in network)x1 structure
    bus         # bus number of wind generator
    Pforecast   # forecasted active power production (pu) 1x(number of time intervals) vector
    Pinj        # actual active power production (pu)
    Pshed       # active power shed from wind forecast through curtailment (pu)

    Wind() = new()
end

type Storage # (number of storage devices in network)x1 structure
    bus             # bus number of storage
    Pdemand         # active power being consumed by device [measured at bus] (pu) 
    powerRate       # power rating (pu)
    energyRate      # energy rating (pu*hr)
    energyLvl       # energy currently stored in device (pu*hr)
    chargeEff       # charging efficiency [between 0 and 1]
    dischargeEff    # discharging efficiency [between 0 and 1]

    Storage() = new()
end

type Params # 1x1 structure
    sBase                   # 3-phase power system base (MVA)
    numOfBus::Integer       # No. of buses in system
    numOfLine::Integer      # No. of lines in system
    numOfGen::Integer       # No. of generators in system
    numOfLoad::Integer      # No. of loads in system
    numOfWind::Integer      # No. of wind generators in system
    numOfStorage::Integer   # No. of storage devices in system
    numOfTimeInt::Integer   # No. of time intervals (for load and wind forecasts)
    hasLosses::Integer      # integer indicating how many segments to use in piecewise-linear model of system losses. (0 indicates that PWL model is not used)
    phaseMax                # maximum phase angle difference allowed across a single line in the network [rad]
    NAI                     # * node arc incidence matrix (nbus x nline, 1 = from -1 = to)
    NGI                     # * node generator incidence matrix (nbus x ngen, 1 = incidence)
    NLI                     # * node load incidence matrix (nbus x nload, 1 = incidence)
    NWI                     # * node wind incidence matrix (nbus x nwind, 1 = incidence)
    NSI                     # * node storage incidence matrix (nbus x nstorage, 1 = incidence)
    phPW_b                  # * PWL approximation of abs(theta_{ij}) for each line

    # Allow partial initialization:
    Params() = new()
end

type Control # 1x1 structure
    T_s             # discrete sampling time (for MPC) [sec] (Must set to 60 for now)
    MPChorizon      # receeding prediction horizon for MPC
    MPCsafeTime     # number of MPC iterations needed without line trip for system to be considered safe 
    MPCmaxTime      # maximum number of MPC iterations allowed
    lineOut         # nline x 1 logical specifying which lines are out of service
    storageTol      # scalar specifying tolerance around 0 for which storage devices are allowed to simulcharge         
    overloadTol     # acceptable line flow overload allowed at safe MPC termination [%]

    Control() = new()
end

type OptParams # 1x1 structure containing GUROBI parameters
    outputflag      # decide whether GUROBI steps should display on console window (1=yes,0=no)
    timelimit       # max time GUROBI is allowed to operate (sec)
    aggregate       # determine aggregator setting (set to 1 to aggregate)

    OptParams() = new()
end

type CostFunction # 1x1 structure containing...
    delThat         # penalty for temperature overload [$/C^2]
    Elvl            # (terminal) cost of deviating from slow-timescale (interpolated) energy storage levels [$/pu^2]
    fG              # (terminal) penalizing generator levels from interpolated references [$/pu^2]
    dfG             # generator ramping cost factor that adds to quadratic gen cost term [$/pu^2]
    windSpill       # penalty for spilling wind from wind farms [$/pu^2]
    loadShed        # penalty for shedding load through demand response [$/pu^2]
    fQ              # storage utility cost [$/pu^2] (i.e. keeping storage level constant is NOT penalized)
    PSTshift        # phase-shift from PST wrt Level 1/2 reference [$/rad^2]

    CostFunction() = new()
end

type Optimize # 1x1 structure
    QPsolver        # optimization solver 2 == GUROBI (1 == CPLEX, currently unavailable)
    params::OptParams
    costFunction::CostFunction

    Optimize() = new()
end

type Time # 1x1 structure
    day     # scalar current day of year
    hour    # scalar current hour of day
    minute  # scalar current minute of hour

    Time() = new()
end

type System
    bus::Bus
    line::Line
    gen::Gen
    load::Load
    wind::Wind
    storage::Storage
    params::Params
    control::Control
    optimize::Optimize
    time::Time

    System() = new()
end

function mat2sys(mpc)
    # This function takes an arbitrary Matpower file as input and creates the
    # basic System structure needed for Level 3.

    baseMVA = mpc["baseMVA"]

    # System.bus
    bus = mpc["bus"]
    
    BUS = Bus(bus[:,1],
            bus[:,2],
            bus[:,10],
            bus[:,8],
            bus[:,9]*π/180.,
            bus[:,5],
            bus[:,6])

    branch = mpc["branch"]
    # System.line
    tap = branch[:,9]
    tap[tap.==0] = 1

    LINE = Line(branch[:,1],
                branch[:,2],
                branch[:,3],
                branch[:,4],
                branch[:,5],
                branch[:,6],
                fill(NaN,length(branch[:,1])),
                zeros(length(branch[:,1])),
                zeros(length(branch[:,1])),
                zeros(length(branch[:,1])),
                zeros(length(branch[:,1])),
                -ones(length(branch[:,1])),
                branch[:,9].!=0,
                branch[:,9].!=0,
                falses(length(branch[:,1])),
                tap,
                branch[:,10].*pi/180,
                zeros(length(branch[:,1])),
                ConductorModel(),
                ThermalModel(),
                fill(NaN,length(branch[:,1])))

    # System.gen
    gen = mpc["gen"]
    Qmin = gen[:,5]
    Qmax = gen[:,4]
    GEN = Gen(
        gen[:,1],
        ones(length(gen[:,1])),
        gen[:,8],
        gen[:,2]./baseMVA,
        gen[:,3]./baseMVA,
        zeros(length(gen[:,1])),
        gen[:,6],
        gen[:,10]./baseMVA,
        gen[:,9]./baseMVA,
        Qmin,#-9999*ones(length(gen[:,1])),
        Qmax,#9999*ones(length(gen[:,1])),
        gen[:,9]./(2*baseMVA),
        -gen[:,9]./(2*baseMVA),
        gen[:,4]./(2*baseMVA),
        gen[:,5]./(2*baseMVA)
        )

    # System.load
    LOAD = Load(
        bus[:,1],
        bus[:,3]./100*ones(1,24),
        bus[:,4]./100*ones(1,24),
        bus[:,3]./baseMVA,
        bus[:,4]./baseMVA,
        zeros(length(bus[:,1])),
        zeros(length(bus[:,1]))
        )

    # System.params
    PARAMS = Params()
    PARAMS.sBase = baseMVA
    PARAMS.numOfBus = length(bus[:,1])
    PARAMS.numOfLine = length(branch[:,1])
    PARAMS.numOfGen = length(gen[:,1])
    PARAMS.numOfLoad = length(bus[:,1])

    # System
    SYSTEM = System()
    SYSTEM.bus = BUS
    SYSTEM.line = LINE
    SYSTEM.gen = GEN
    SYSTEM.load = LOAD
    SYSTEM.params = PARAMS

    return SYSTEM
end

end