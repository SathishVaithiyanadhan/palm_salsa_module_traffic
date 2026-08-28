{salsa+simple_tra.spc
Current revision
----------------
 20180903 Added SALSA variables               monakurppa
 20260826 Created traffic tagging, H2SO4 equation, OCSV formation  Sathish
Former revisions
----------------
 $Id: smog.spc 2459 2017-09-13 14:10:33Z forkel $
}
#include atoms

 #DEFVAR
    O       = O ;      	    {oxygen atomic ground state (3P)}
    O3	    = 3O ;          {ozone}  
    NO	    = N + O ;       {nitric oxide}  
    NO2	    = N + 2O ;      {nitrogen dioxide} 
    NO3	    = N + 3O ;      {nitrogen trioxide} 
    N2O5    = 2N + 5O ;     {dinitrogen pentoxide}           
    HNO3    = H + N + 3O ;  { nitric acid }
    HNO4    = H + N + 4O ;  {HO2NO2 pernitric acid}
    H       = H ;           {hydrogen atomic ground state (2S)}
    OH	    = O + H ;       {hydroxyl radical}  
    HO2	    = H + 2O ;      {perhydroxyl radical}                     
    H2O2    = 2H + 2O ;     {hydrogen peroxide} 
    CH3     = C + 3H ;      {methyl radical}
    CH3O    = C + 3H + O ;  {methoxy radical}
    CH3O2   = C + 3H + 2O ; {methylperoxy radical}
    CH3OOH  = C + 4H + 2O ; {CH4O2 methylperoxy alcohol}
    HCO     = H + C + O ;   {CHO formyl radical}
    CH2O    = C + 2H + O ;  {formalydehyde}
    CO      = C + O ;       {carbon monoxide}   
    
    RH      = ignore ;      {alkanes} 
    RO2     = ignore ;      {alkyl peroxy radical}
    RCHO    = ignore ;      {carbonyl}
    RCOO2   = ignore ;      
    RCOO2NO2= ignore ;

    H2SO4   = 2H + S +4O ;  {sulfuric acid}
    SO2     = S + 2O ;        {sulfur dioxide}
    NH3     = 3H + N ;      {ammonia}
    OCNV    = ignore ;      {non-volatile OC}
    OCSV    = ignore ;      {semi-volatile OC}
    
      // =====  TRAFFIC-TAGGED SPECIES =====
      // NOx tagging with Ox/NOy split (TOAST 2018):
      //   NO2_tra  (NOy) = null-cycle reservoir, never creates traffic O3
      //   NO2_x_tra(Ox)  = net-production conduit to traffic O3

    O3_tra   = 3O ;         {traffic-tagged ozone (Ox)}
    NO_tra   = N + O ;      {traffic-tagged nitric oxide (NOy)}  
    NO2_tra  = N + 2O ;     {traffic-tagged NO2 (NOy, null-cycle reservoir)}
    NO2_x_tra = N + 2O ;    {traffic-tagged NO2 (Ox, net production -> O3_tra)}
    OH_tra   = O + H ;      {traffic-tagged OH (diagnostic tracer)}
    HO2_tra  = H + 2O ;     {traffic-tagged perhydroxyl radical (chain)}
    RH_tra   = ignore ;     {traffic-tagged alkanes}
    RO2_tra  = ignore ;     {traffic-tagged alkyl peroxy radical}
    RCHO_tra = ignore ;     {traffic-tagged carbonyl}
    HNO3_tra = H + N + 3O ; {traffic-tagged nitric acid}
    H2SO4_tra = 2H + S + 4O ; {traffic-tagged sulfuric acid}
    SO2_tra     = S + 2O ;        {traffic-tagged sulfur dioxide}
    NH3_tra  = 3H + N ;     {traffic-tagged ammonia}
    OCNV_tra = ignore ;     {traffic-tagged non-volatile OC}
    OCSV_tra = ignore ;     {traffic-tagged semi-volatile OC}
    
   // =====  PASSIVE PARTICULATE MATTER =====
    PM10     = ignore ;     {PM10 - particulate matter < 10 microns, passive tracer}
    PM25     = ignore ;     {PM2.5 - particulate matter < 2.5 microns, passive tracer}
    PM10_tra = ignore ;     {traffic-tagged PM10 - passive tracer}
    PM25_tra = ignore ;     {traffic-tagged PM2.5 - passive tracer}

#DEFFIX
    H2O     = H + 2O ;      {water}
    H2      = 2H ;          {molecular hydrogen}
    O2      = 2O ;          {molecular oxygen}
    N2      = 2N ;          {molecular nitrogen}
    CH4     = C + 4H ;      {methane}
    CO2     = C + 2O ;      {carbon dioxide}