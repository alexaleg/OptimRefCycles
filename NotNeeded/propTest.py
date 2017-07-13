import CoolProp as CP

HEOS = CP.AbstractState("HEOS", "Methane&Ethane")

HEOS.set_mole_fractions([0.2, 0.8])
HEOS.update(CP.PT_INPUTS, 101325, 200)
HEOS.mole_fractions_liquid()
HEOS.mole_fractions_vapor()
HEOS.Q()
