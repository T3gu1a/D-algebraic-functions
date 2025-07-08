#Arithmetic of Difference-Algebraic Sequences

DalgSeq := module()

option `Copyright (c) 2024 Bertrand Teguia Tabuguia, University of Oxford`, package;

export arithmeticDalgSeq, unaryDalgSeq, PartialSumDalgSeq, PartialProdDalgSeq, 
       radicalDalgSeq, CCfiniteToDalg, HoloToSimpleRatrec, OrderDegreeRec, DDStoADE,
       DalgGuess, deltakshift, AnsatzDalgSeq;

local buildsystemseq, mergesystemseq, SystoDE, HoloToSimpleRatrecLA, HoloToSimpleRatrecGB,
      checkingguess, REorders, shiftkstart;

$include <NLDE/DalgSeq/arithmeticDalgSeq/src/arithmeticDalgSeq.mm>
$include <NLDE/DalgSeq/unaryDalgSeq/src/unaryDalgSeq.mm>
$include <NLDE/DalgSeq/radicalDalgSeq/src/radicalDalgSeq.mm>
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
$include <NLDE/DalgSeq/DalgGuess/src/DalgGuess.mm>
$include <NLDE/DalgSeq/DalgGuess/SubProcedures/src/deltakshift.mm>
$include <NLDE/DalgSeq/DalgGuess/SubProcedures/src/checkingguess.mm>
$include <NLDE/DalgSeq/DalgGuess/SubProcedures/src/shiftkstart.mm>
$include <NLDE/DalgSeq/AnsatzDalgSeq/src/AnsatzDalgSeq.mm>
$include <NLDE/DalgSeq/CCfiniteToSimpleRatrec/src/CCfiniteToSimpleRatrec.mm>

end module: #end DalgSeq