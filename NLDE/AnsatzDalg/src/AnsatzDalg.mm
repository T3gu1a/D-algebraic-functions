
AnsatzDalg:= module()

option `Copyright (c) 2022 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

export deltakdiff, unaryDeltak, arithmeticDeltak;

local buildsystem, mergesystem, ComputDegkDE, DegreekDE, startkorder, ordertoktuple;

$include <NLDE/AnsatzDalg/unaryDeltak/src/unaryDeltak.mm>
$include <NLDE/AnsatzDalg/arithmeticDeltak/src/arithmeticDeltak.mm>
$include <NLDE/AnsatzDalg/CommonInternalProcedures/src/mergesystem.mm>
$include <NLDE/AnsatzDalg/CommonInternalProcedures/src/DegreekDE.mm>

end module: #end AnsatzDalg