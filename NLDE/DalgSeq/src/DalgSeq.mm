#Arithmetic of Difference-Algebraic Sequences

DalgSeq := module()

option `Copyright (c) 2024 Bertrand Teguia Tabuguia, University of Oxford`, package;

export arithmeticDalgSeq, unaryDalgSeq, FiniteSeriesDalgSeq, FiniteProdDalgSeq, 
       REorders, CCfiniteToDalg, HoloToSimpleRatrec, OrderDegreeRec, DDStoADE;

local buildsystemseq, mergesystemseq, SystoDE, HoloToSimpleRatrecLA, HoloToSimpleRatrecGB;

$include <NLDE/DalgSeq/arithmeticDalgSeq/src/arithmeticDalgSeq.mm>
$include <NLDE/DalgSeq/unaryDalgSeq/src/unaryDalgSeq.mm>
$include <NLDE/DalgSeq/FiniteSeriesDalgSeq/src/FiniteSeriesDalgSeq.mm>
$include <NLDE/DalgSeq/FiniteProdDalgSeq/src/FiniteProdDalgSeq.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/buildsystemseq.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/mergesystemseq.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/SystoDE.mm>
$include <NLDE/DalgSeq/CommonInternalProcedures/src/REorders.mm>
$include <NLDE/DalgSeq/CCfiniteToDalg/src/CCfiniteToDalg.mm>
$include <NLDE/DalgSeq/HoloToSimpleRatrec/src/HoloToSimpleRatrec.mm>
$include <NLDE/DalgSeq/OrderDegreeRec/src/OrderDegreeRec.mm>
$include <NLDE/DalgSeq/DDStoADE/src/DDStoADE.mm>

end module: #end DalgSeq