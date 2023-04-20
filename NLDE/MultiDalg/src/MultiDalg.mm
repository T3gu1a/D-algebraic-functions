MultiDalg := module()

option `Copyright (c) 2023 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

export thetaderiv, arithmeticMDalg, CantorSigma, CantorInvSigma;

local   deglexisless, thetakTuple::table(table):=table([]), sortedkpartition, compwiseless, derivation;

$include <NLDE/MultiDalg/arithmeticMDalg/src/arithmeticMDalg.mm>
$include <NLDE/MultiDalg/thetaderiv/src/thetaderiv.mm>
$include <NLDE/MultiDalg/CantorSigma/src/CantorSigma.mm>
$include <NLDE/MultiDalg/CantorInvSigma/src/CantorInvSigma.mm>
$include <NLDE/MultiDalg/CommonInternalProcedures/src/deglexisless.mm>
$include <NLDE/MultiDalg/CommonInternalProcedures/src/sortedkpartition.mm>
$include <NLDE/MultiDalg/CommonInternalProcedures/src/derivation.mm>
$include <NLDE/MultiDalg/CommonInternalProcedures/src/compwiseless.mm>	

end module: #end MultiDalg