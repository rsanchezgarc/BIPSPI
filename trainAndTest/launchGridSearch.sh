#sampledComplexesPath="~/Tesis/rriPredMethod/data/bench5Data/newCodeData/codifiedInput/struct_2/sampledInputs"
sampledComplexesPath="~/Tesis/rriPredMethod/data/ppdockingBenchData/newCodeData/codifiedInput/struct_2/sampledInputs"
numProc=2
nFolds=3
python -m trainAndTest.trainAndTest -g -v -j $numProc -k $nFolds -i $sampledComplexesPath
