
       module pac_data

       implicit none

c      -- ELEMENTS DATA
c       - Atomic weights
c      ...nmaxelem=72     in dimension.common !!!
       character*2, public, dimension (72)  :: element        ! element name
       real*8, public, dimension (72)       :: atomic_weight  ! atomic weight [amu]

c      -- ATOMS AND MOLECULES DATA
c       - Lennard-Jones parameters to estimate the coefficient of molecular diffusion D
c      ...nmaxlj=52       in dimension.common !!!
c      ...nmaxcharspec=12 in dimension.common !!!
       character*12, public, dimension (52) :: lj_spec        ! species name
       real*8, public, dimension (52)       :: lj_length      ! Lennard-Jones characteristic length [A]
       real*8, public, dimension (52)       :: lj_energy      ! Lennard-Jones characteristic energy [K]
cc       - Polarizability to estimate the Rayleigh scattering cross section (read in species file !)
cc      ...nmaxpo=96       in dimension.common !!!
cc      ...nmaxcharspec=12 in dimension.common !!!
c       character*12, public, dimension (96) :: po_spec        ! species name
c       real*8, public, dimension (96)       :: polary         ! polarizability [A^3]

       data element/
     #    'H ',
     #    'He',
     #    'Li',
     #    'Be',
     #    'B ',
     #    'C ',
     #    'N ',
     #    'O ',
     #    'F ',
     #    'Ne',
     #    'Na',
     #    'Mg',
     #    'Al',
     #    'Si',
     #    'P ',
     #    'S ',
     #    'Cl',
     #    'Ar',
     #    'K ',
     #    'Ca',
     #    'Sc',
     #    'Ti',
     #    'V ',
     #    'Cr',
     #    'Mn',
     #    'Fe',
     #    'Co',
     #    'Ni',
     #    'Cu',
     #    'Zn',
     #    'Ga',
     #    'Ge',
     #    'As',
     #    'Se',
     #    'Br',
     #    'Kr',
     #    'Rb',
     #    'Sr',
     #    'Y ',
     #    'Zr',
     #    'Nb',
     #    'Mo',
     #    'Tc',
     #    'Ru',
     #    'Rh',
     #    'Pd',
     #    'Ag',
     #    'Cd',
     #    'In',
     #    'Sn',
     #    'Sb',
     #    'Te',
     #    'I ',
     #    'Xe',
     #    'Cs',
     #    'Ba',
     #    'Lu',
     #    'Hf',
     #    'Ta',
     #    'W ',
     #    'Re',
     #    'Os',
     #    'Ir',
     #    'Pt',
     #    'Au',
     #    'Hg',
     #    'Tl',
     #    'Pb',
     #    'Bi',
     #    'Po',
     #    'At',
     #    'Rn'/
       data atomic_weight/ 
     #    1.00794d0,    ! H  Z=1
     #    4.002602d0,   ! He Z=2
     #    6.941d0,      ! Li Z=3
     #    9.012182d0,   ! Be Z=4
     #   10.811d0,      ! B  Z=5
     #   12.0107d0,     ! C  Z=6
     #   14.0067d0,     ! N  Z=7
     #   15.9994d0,     ! O  Z=8
     #   18.9984032d0,  ! F  Z=9
     #   20.1797d0,     ! Ne Z=10
     #   22.98976928d0, ! Na Z=11
     #   24.3050d0,     ! Mg Z=12
     #   26.9815386d0,  ! Al Z=13
     #   28.0855d0,     ! Si Z=14
     #   30.973762d0,   ! P  Z=15
     #   32.065d0,      ! S  Z=16
     #   35.453d0,      ! Cl Z=17
     #   39.948d0,      ! Ar Z=18
     #   39.0983d0,     ! K  Z=19
     #   40.078d0,      ! Ca Z=20
     #   44.955912d0,   ! Sc Z=21
     #   47.867d0,      ! Ti Z=22
     #   50.9415d0,     ! V  Z=23
     #   51.9961d0,     ! Cr Z=24
     #   54.938045d0,   ! Mn Z=25
     #   55.845d0,      ! Fe Z=26
     #   58.933195d0,   ! Co Z=27
     #   58.6934d0,     ! Ni Z=28
     #   63.546d0,      ! Cu Z=29
     #   65.38d0,       ! Zn Z=30
     #   69.723d0,      ! Ga Z=31
     #   72.64d0,       ! Ge Z=32
     #   74.92160d0,    ! As Z=33
     #   78.96d0,       ! Se Z=34
     #   79.904d0,      ! Br Z=35
     #   83.798d0,      ! Kr Z=36
     #   85.4678d0,     ! Rb Z=37
     #   87.62d0,       ! Sr Z=38
     #   88.90585d0,    ! Y  Z=39
     #   91.224d0,      ! Zr Z=40
     #   92.90638d0,    ! Nb Z=41
     #   95.96d0,       ! Mo Z=42
     #   98.d0,         ! Tc Z=43
     #   101.07d0,      ! Ru Z=44
     #   102.90550d0,   ! Rh Z=45
     #   106.42d0,      ! Pd Z=46
     #   107.8682d0,    ! Ag Z=47
     #   112.411d0,     ! Cd Z=48
     #   114.818d0,     ! In Z=49
     #   118.710d0,     ! Sn Z=50
     #   121.760d0,     ! Sb Z=51
     #   127.60d0,      ! Te Z=52
     #   126.90447d0,   ! I  Z=53
     #   131.293d0,     ! Xe Z=54
     #   132.9054519d0, ! Cs Z=55
     #   137.327d0,     ! Ba Z=56
     #   174.9668d0,    ! Lu Z=71
     #   178.49d0,      ! Hf Z=72
     #   180.94788d0,   ! Ta Z=73
     #   183.84d0,      ! W  Z=74
     #   186.207d0,     ! Re Z=75
     #   190.23d0,      ! Os Z=76
     #   192.217d0,     ! Ir Z=77
     #   195.084d0,     ! Pt Z=78
     #   196.966569d0,  ! Au Z=79
     #   200.59d0,      ! Hg Z=80
     #   204.3833d0,    ! Tl Z=81
     #   207.2d0,       ! Pb Z=82
     #   208.98040d0,   ! Bi Z=83
     #   209.d0,        ! Po Z=84
     #   210.d0,        ! At Z=85
     #   222.d0/        ! At Z=86

       data lj_spec/
c      ...values from Reid et al. 1980, "The Properties of Gases and Liquids", McGraw Hill, Appendix B
     # 'Ar          ','He          ','Kr          ','Ne          ',
     # 'CH3OH       ','CH4         ','CO          ','OCS         ',
     # 'CO2         ','CS2         ','C2H2        ','C2H4        ',
     # 'C2H6        ','C2H5OH      ','NCCN        ','CH3OCH3     ',
     # 'C3H6        ','p-C3H4      ','C3H8        ','CH3COCH3    ',
     # 'n-C4H10     ','i-C4H10     ','b-C6H6      ','Cl2         ',
     # 'F2          ','HBr         ','HCN         ','HCl         ',
     # 'HF          ','HI          ','H2          ','H2O         ',
     # 'H2O2        ','H2S         ','I2          ','NH3         ',
     # 'NO          ','N2          ','N2O         ','O2          ',
     # 'PH3         ','SO2         ','SiF4        ','SiH4        ',
c      ...lj_length values for atoms: van der Waals radius to the power of 1.65 (estimated from H2,N2,O2,Cl2)
     # 'H           ','C           ','N           ','O           ',
     # 'F           ','P           ','S           ','Cl          '/
       data lj_length/
     # 3.542d0       ,2.551d0       ,3.655d0       ,2.820d0       ,
     # 3.626d0       ,3.758d0       ,3.69d0        ,4.13d0        ,
     # 3.941d0       ,4.483d0       ,4.033d0       ,4.163d0       ,
     # 4.443d0       ,4.530d0       ,4.361d0       ,4.307d0       ,
     # 4.678d0       ,4.761d0       ,5.118d0       ,4.600d0       ,
     # 4.687d0       ,5.278d0       ,5.349d0       ,4.217d0       ,
     # 3.357d0       ,3.353d0       ,3.630d0       ,3.339d0       ,
     # 3.148d0       ,4.211d0       ,2.827d0       ,2.641d0       ,
     # 4.196d0       ,3.623d0       ,5.160d0       ,2.900d0       ,
     # 3.492d0       ,3.798d0       ,3.828d0       ,3.467d0       ,
     # 3.981d0       ,4.112d0       ,4.880d0       ,4.084d0       ,
c
     # 1.35d0        ,2.40d0        ,2.06d0        ,2.00d0        ,
     # 1.89d0        ,2.64d0        ,2.64d0        ,2.52d0        /
       data lj_energy/
     # 93.3d0        ,10.22d0       ,178.9d0       ,32.8d0        ,
     # 481.4d0       ,148.6d0       ,91.7d0        ,336.0d0       ,
     # 195.2d0       ,467.0d0       ,231.8d0       ,224.7d0       ,
     # 215.7d0       ,362.6d0       ,348.6d0       ,395.0d0       ,
     # 248.9d0       ,251.8d0       ,237.1d0       ,560.2d0       ,
     # 531.4d0       ,330.1d0       ,412.3d0       ,316.0d0       ,
     # 112.6d0       ,449.0d0       ,569.1d0       ,344.7d0       ,
     # 330.0d0       ,288.7d0       ,59.7d0        ,809.1d0       ,
     # 289.3d0       ,301.1d0       ,474.2d0       ,558.3d0       ,
     # 116.7d0       ,71.4d0        ,232.4d0       ,106.7d0       ,
     # 251.5d0       ,335.4d0       ,171.9d0       ,207.6d0       ,
c
     # -1.0d0        ,-1.0d0        ,-1.0d0        ,-1.0d0        ,
     # -1.0d0        ,-1.0d0        ,-1.0d0        ,-1.0d0        /

c       data po_spec/
cc      ...polarizability values from: Lide, R D, CRC Handbook of Chemistry and Physics, 90th ed, CRC Press 2009 (experimental)
cc                                     Woon & Herbst 2009, ApJS, 185, 273 (computational)
c     # 'H           ','H2          ','He          ','C           ',
c     # 'O           ','N           ','S           ','Si          ',
c     # 'P           ','F           ','Cl          ','Fe          ',
c     # 'Na          ','Mg          ','CH          ','CH2         ',
c     # 'CH3         ','CH4         ','C2          ','C2H         ',
c     # 'C2H2        ','C2H3        ','C2H4        ','C2H5        ',
c     # 'C2H6        ','C3          ','C-C3H       ','C3H         ',
c     # 'C-C3H2      ','H2C3        ','C3H2        ','C3H3        ',
c     # 'P-C3H4      ','A-C3H4      ','C4          ','C4H         ',
c     # 'C4H2        ','N-C4H3      ','I-C4H3      ','OH          ',
c     # 'H2O         ','O2          ','HO2         ','H2O2        ',
c     # 'CO          ','HCO         ','H2CO        ','CH2OH       ',
c     # 'CH3O        ','CH3OH       ','CO2         ','HCOOH       ',
c     # 'C2O         ','CH2CO       ','CH3CHO      ','C3O         ',
c     # 'NH          ','NH2         ','NH3         ','N2          ',
c     # 'CN          ','HCN         ','HNC         ','H2CN        ',
c     # 'CH2NH       ','C2N         ','C3N         ','HC3N        ',
c     # 'CH2CN       ','CH3CN       ','C2H3CN      ','NH2CN       ',
c     # 'NCCN        ','CH3C3N      ','NO          ','HNO         ',
c     # 'NCO         ','N2O         ','NO2         ','SH          ',
c     # 'H2S         ','S2          ','HS2         ','H2S2        ',
c     # 'CS          ','HCS         ','H2CS        ','CS2         ',
c     # 'C2S         ','C3S         ','SO          ','SO2         ',
c     # 'OCS         ','NS          ','HF          ','HCl         '/
c       data polary/
c     # 0.666793d0    ,0.804d0       ,0.2050522d0   ,1.76d0        , ! Lide 2009 (H2: Woon & Herbst 2009)
c     # 0.802d0       ,1.10d0        ,2.90d0        ,5.38d0        , ! Lide 2009
c     # 3.63d0        ,0.557d0       ,2.18d0        ,8.4d0         , ! Lide 2009
c     # 24.11d0       ,10.6d0        ,2.12d0        ,2.12d0        , ! Lide 2009 (CH,CH: Woon & Herbst 2009)
c     # 2.335d0       ,2.593d0       ,5.074d0       ,4.415d0       , ! Woon & Herbst 2009 (CH4: Lide 2009)
c     # 3.487d0       ,3.856d0       ,4.188d0       ,4.103d0       , ! Woon & Herbst 2009
c     # 4.47d0        ,5.179d0       ,4.802d0       ,5.359d0       , ! Woon & Herbst 2009 (C2H6: Lide 2009)
c     # 4.583d0       ,5.609d0       ,2.581d0       ,5.556d0       , ! Woon & Herbst 2009 (C3H3: The Astrochymist web)
c     # 6.18d0        ,6.18d0        ,7.512d0       ,7.151d0       , ! Woon & Herbst 2009 (A-C3H4,P-C3H4: Lide 2009)
c     # 7.048d0       ,7.643d0       ,7.643d0       ,1.073d0       , ! Woon & Herbst 2009 (N-C4H3,I-C4H3: c3ch3 from The Astrochymist web)
c     # 1.429d0       ,1.569d0       ,1.954d0       ,2.243d0       , ! Woon & Herbst 2009
c     # 1.953d0       ,2.538d0       ,2.770d0       ,3.0d0         , ! Woon & Herbst 2009 (CH2OH: assumed)
c     # 3.0d0         ,3.32d0        ,2.665d0       ,3.317d0       , ! Woon & Herbst 2009 (CH3O: assumed, CH3OH: Lide 2009)
c     # 4.087d0       ,4.301d0       ,4.278d0       ,6.027d0       , ! Woon & Herbst 2009
c     # 1.447d0       ,1.782d0       ,2.22d0        ,1.74d0        , ! Woon & Herbst 2009 (NH3: Lide 2009)
c     # 2.884d0       ,2.53d0        ,2.729d0       ,2.963d0       , ! Woon & Herbst 2009
c     # 3.327d0       ,4.270d0       ,5.675d0       ,5.848d0       , ! Woon & Herbst 2009
c     # 4.382d0       ,4.48d0        ,6.109d0       ,3.947d0       , ! Woon & Herbst 2009
c     # 7.99d0        ,8.008d0       ,1.698d0       ,2.151d0       , ! Woon & Herbst 2009 (NCCN: Lide 2009)
c     # 3.176d0       ,2.919d0       ,2.910d0       ,3.261d0       , ! Woon & Herbst 2009
c     # 3.631d0       ,5.994d0       ,6.474d0       ,6.567d0       , ! Woon & Herbst 2009
c     # 4.262d0       ,4.882d0       ,5.014d0       ,8.74d0        , ! Woon & Herbst 2009 (CS2: Lide 2009)
c     # 6.873d0       ,9.649d0       ,3.422d0       ,3.219d0       , ! Woon & Herbst 2009
c     # 4.67d0        ,3.785d0       ,0.80d0        ,2.515d0       / ! Woon & Herbst 2009c

       end module pac_data
