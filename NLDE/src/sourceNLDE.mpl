NLDE:= module()

option `Copyright (c) 2022 Bertrand Teguia Tabuguia, Max Planck Institute for MiS, Leipzig`, package;

export unaryDalg, diffDalg, invDalg, SysToMinDiffPoly, composeDalg, arithmeticDalg, AnsatzDalg, DDfiniteToDalg;

local buildsystem, mergesystem, ftogh, subsgfurther, ftogx, NLDE_nlho;

$include <NLDE/CommonInternalProcedures/src/mergesystem.mm>
$include <NLDE/SysToMinDiffPoly/src/SysToMinDiffPoly.mm>
$include <NLDE/composeDalg/src/composeDalg.mm>
$include <NLDE/unaryDalg/src/unaryDalg.mm>	
$include <NLDE/arithmeticDalg/src/arithmeticDalg.mm>
$include <NLDE/diffDalg/src/diffDalg.mm>	
$include <NLDE/invDalg/src/invDalg.mm>
$include <NLDE/DDfiniteToDalg/src/DDfiniteToDalg.mm>
$include <NLDE/CommonInternalProcedures/src/NLDE_nlho.mm>
$include <NLDE/AnsatzDalg/src/AnsatzDalg.mm>

end module:

savelib('NLDE',"path_to_your_folder_for_softwarepackages/NLDE.mla"):
