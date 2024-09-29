#Arithmetic of Difference-Algebraic Sequences

DalgSeq := module()

option `Copyright (c) 2024 Bertrand Teguia Tabuguia, University of Oxford`, package;

export arithmeticDalgSeq, unaryDalgSeq, FiniteSeriesDalgSeq, FiniteProdDalgSeq, REorders;

local buildsystemseq, mergesystemseq, SystoDE;

$include <NLDE/DalgSeq/arithmeticDalgSeq/src/arithmeticDalgSeq.mm>
$include <NLDE/DalgSeq/unaryDalgSeq/src/unaryDalgSeq.mm>
$include <NLDE/DalgSeq/FiniteSeriesDalgSeq/src/FiniteSeriesDalgSeq.mm>
$include <NLDE/DalgSeq/FiniteProdDalgSeq/src/FiniteProdDalgSeq.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/buildsystemseq.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/mergesystemseq.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/SystoDE.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/REorders.mm>

end module: #end DalgSeq