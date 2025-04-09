
AnsatzDalgSeq:= module()

option `Copyright (c) 2025 Bertrand Teguia Tabuguia, University of Oxford`, package;

export unaryDeltakSeq, arithmeticDeltakSeq;

local buildsystem, mergesystem, ComputDegkDE, DegreekDE;

$include <NLDE/DalgSeq/AnsatzDalgSeq/unaryDeltakSeq/src/unaryDeltakSeq.mm>
$include <NLDE/DalgSeq/AnsatzDalgSeq/arithmeticDeltakSeq/src/arithmeticDeltakSeq.mm>
$include <NLDE/DalgSeq/AnsatzDalgSeq/CommonInternalProcedures/src/mergesystem.mm>
$include <NLDE/DalgSeq/AnsatzDalgSeq/CommonInternalProcedures/src/DegreekDE.mm>

end module: #end AnsatzDalgSeq