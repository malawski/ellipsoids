var argv = require('optimist').argv;

// this signals creates prerequisites
function createBasicSignals(wfOut, pathToRootDir, pathToEllipsoidProgram) {
	var signals = []
	signals.push("rootPath")
	signals.push("ellipsoidPath")
	signals.push("datDirs")
	signals.push("endMessage")

	for (i = 0; i < signals.length; i++) {
		wfOut.signals.push({
			name : signals[i]
		});
	}

	wfOut.signals[0].data = [ pathToRootDir ]
	wfOut.signals[1].data = [ pathToEllipsoidProgram ]
	wfOut.signals[2].data = []
	wfOut.signals[3].data = []
}
// this funciton create signals, as they can be passed between processes with two indexes
function createUsualSignals(howMuch, indexOfBaseExecution, prefix_name, wfOut) {

	for (i = 0; i < howMuch; i++) {
		wfOut.signals.push({
			name : prefix_name + i + "_" + indexOfBaseExecution,
			data : []
		});
	}
}

//this funciton create signals, as they can be passed between processes, with one index
function createUsualSimpleSignals(howMuch, prefix_name, wfOut) {

	for (i = 0; i < howMuch; i++) {
		wfOut.signals.push({
			name : prefix_name + i,
			data : []
		});
	}
}

// this function creates data for functionExecution processes
function createSignalDat(wfOut, index, ifCoord, numberOfParticles,
		numberOfSpecies, forceScalingFactor, rotationScalingFactor,
		diameterIncreasingFactor, cellsX, cellsY, cellsZ, constractionRate,
		initialPackingDensity, maximalNumberOfSteps, isConstantVolume,
		numberOfParts, diameterOfParts1, diameterOfParts2, diameterOfParts3,
		numberOfLinesPerPage, numberOfStepsBetweenPrintouts,
		numberOfStepsBetweenCoord, numberOfStepsBetweenRotations) {

	return {
		"name" : "tobechanged",
		"ifCoord" : ifCoord,
		"numberOfParticles" : numberOfParticles,
		"numberOfSpecies" : numberOfSpecies,
		"forceScalingFactor" : forceScalingFactor,
		"rotationScalingFactor" : rotationScalingFactor,
		"diameterIncreasingFactor" : diameterIncreasingFactor,
		"cellsX" : cellsX,
		"cellsY" : cellsY,
		"cellsZ" : cellsZ,
		"constractionRate" : constractionRate,
		"initialPackingDensity" : initialPackingDensity,
		"maximalNumberOfSteps" : maximalNumberOfSteps,
		"isConstantVolume" : isConstantVolume,
		"numberOfParts" : numberOfParts,
		"diameterOfParts1" : diameterOfParts1,
		"diameterOfParts2" : diameterOfParts2,
		"diameterOfParts3" : diameterOfParts3,
		"numberOfLinesPerPage" : numberOfLinesPerPage,
		"numberOfStepsBetweenPrintouts" : numberOfStepsBetweenPrintouts,
		"numberOfStepsBetweenCoord" : numberOfStepsBetweenCoord,
		"numberOfStepsBetweenRotations" : numberOfStepsBetweenRotations
	}

}

// this function creates prcess which will execute console program
function task(name, functionName, executable, args, ins, outs) {
	return {
		"name" : name,
		"function" : functionName,
		"type" : "dataflow",
		// "firingLimit": 1,
		"config" : {
			"executor" : {
				// "queue_name": "test1",
				"executable" : executable,
				"args" : args
			}
		},
		"ins" : ins,
		"outs" : outs
	}
}

// this function creates simple process which refers to function in functions.js
// file
function functionExecution(name, functionName, ins, outs) {
	return {
		"name" : name,
		"function" : functionName,
		"type" : "dataflow",
		"ins" : ins,
		"outs" : outs
	}
}
function createWf(functionName, pathToRootDir, pathToEllipsoidProgram,
		numberOfBaseExecutions, pathToAveragingProgram) {
	// main tuple
	var wfOut = {
		processes : [],
		signals : [],
		ins : [],
		outs : [ "endMessage" ]
	};

	// data variables
	var ifCoord = [ "1" ];
	var numberOfParticles = [ "1000" ];
	var numberOfSpecies = [ "1" ];
	var forceScalingFactor = [ "0.1" ];
	var rotationScalingFactor = [ "3.00" ];
	var diameterIncreasingFactor = [ "0.01" ];
	var cellsX = [ "20", "40" ];
	var cellsY = [ "20" ];
	var cellsZ = [ "20" ];
	var constractionRate = [ "56" ];
	var initialPackingDensity = [ "0.8" ];
	var maximalNumberOfSteps = [ "10000000" ];
	var isConstantVolume = [ "1" ];
	var numberOfParts = [ "1000" ];
	var diameterOfParts1 = [ "1.3" ];
	var diameterOfParts2 = [ "1.0" ];
	var diameterOfParts3 = [ "1.0" ];
	var numberOfLinesPerPage = [ "56" ];
	var numberOfStepsBetweenPrintouts = [ "100" ];
	var numberOfStepsBetweenCoord = [ "1000000" ];
	var numberOfStepsBetweenRotations = [ "1" ];

	var lengths = [ ifCoord.length, numberOfParticles.length,
			numberOfSpecies.length, forceScalingFactor.length,
			rotationScalingFactor.length, diameterIncreasingFactor.length,
			cellsX.length, cellsY.length, cellsZ.length,
			constractionRate.length, initialPackingDensity.length,
			maximalNumberOfSteps.length, isConstantVolume.length,
			numberOfParts.length, diameterOfParts1.length,
			diameterOfParts2.length, diameterOfParts3.length,
			numberOfLinesPerPage.length, numberOfStepsBetweenPrintouts.length,
			numberOfStepsBetweenCoord.length,
			numberOfStepsBetweenRotations.length ];

	var counters = [ lengths.length ];
	var product = 1;
	for (i = 0; i < lengths.length; i++) {
		counters[i] = lengths[i];
		product *= lengths[i];
	}
	var pathToRoot = String(pathToRootDir)
	var pathToAveraging = String(pathToAveragingProgram)
	var iterationNumber = 0
	// ins
	wfOut.ins.push("rootPath");
	
	// signals

	createBasicSignals(wfOut, pathToRootDir, pathToEllipsoidProgram);
	// loop below creates signal with data for all permutation
	for (i = 0; i < product; i++) {
		var datDirName = pathToRoot + '/' + i;
		// console.log("This is path to dat dir name " + datDirName);
		var repeatNumber = 0;
		for (repeatNumber = 0; repeatNumber < numberOfBaseExecutions; repeatNumber++) {
			var index = 0;
			var mySignal = [];
			mySignal = createSignalDat(wfOut, iterationNumber,
					ifCoord[counters[index++] - 1],
					numberOfParticles[counters[index++] - 1],
					numberOfSpecies[counters[index++] - 1],
					forceScalingFactor[counters[index++] - 1],
					rotationScalingFactor[counters[index++] - 1],
					diameterIncreasingFactor[counters[index++] - 1],
					cellsX[counters[index++] - 1],
					cellsY[counters[index++] - 1],
					cellsZ[counters[index++] - 1],
					constractionRate[counters[index++] - 1],
					initialPackingDensity[counters[index++] - 1],
					maximalNumberOfSteps[counters[index++] - 1],
					isConstantVolume[counters[index++] - 1],
					numberOfParts[counters[index++] - 1],
					diameterOfParts1[counters[index++] - 1],
					diameterOfParts2[counters[index++] - 1],
					diameterOfParts3[counters[index++] - 1],
					numberOfLinesPerPage[counters[index++] - 1],
					numberOfStepsBetweenPrintouts[counters[index++] - 1],
					numberOfStepsBetweenCoord[counters[index++] - 1],
					numberOfStepsBetweenRotations[counters[index++] - 1]);

			mySignal.name = "signal_" + iterationNumber + "_" + repeatNumber;
			wfOut.signals.push(mySignal);

		}
		// permutation index
		iterationNumber += 1

		counters[lengths.length - 1]--;
		for (j = lengths.length - 1; j > 0; j--) {
			if (counters[j] == 0) {
				counters[j] = lengths[j];
				counters[j - 1]--;
			}
		}
	}
	for (j = 0; j < numberOfBaseExecutions; j++) {
		createUsualSignals(iterationNumber, j, "path_", wfOut);
	}
	for (j = 0; j < numberOfBaseExecutions; j++) {
		createUsualSignals(iterationNumber, j, "done_", wfOut);
	}
	var signalNames = []
	for (i = 0; i < iterationNumber; i++) {
		for (j = 0; j < numberOfBaseExecutions; j++) {
			signalNames.push("signal_" + i + "_" + j)
		}
	}

	createUsualSimpleSignals(iterationNumber, "average_done_", wfOut);
	createUsualSimpleSignals(iterationNumber, "signal_permutation_average_", wfOut);
	// processes
	wfOut.processes.push(functionExecution("mkdir", "createRootDir", [ 0,
			"rootPath" ], signalNames)); // this makes root directory

	for (i = 0; i < iterationNumber; i++) {
		for (j = 0; j < numberOfBaseExecutions; j++) {
			wfOut.processes.push(functionExecution("create_dat_" + i + "_" + j,
					"createDat", [ "signal_" + i + "_" + j, "rootPath" ], [
							"path_" + i + "_" + j, "datDirs" ])); // this
																	// creates
																	// local
																	// directories
																	// and
		} // fulfils directories with dat files
	}
	for (i = 0; i < iterationNumber; i++) {
		for (j = 0; j < numberOfBaseExecutions; j++) {
			var pathToDatFile = pathToRootDir + "/" + i + "/" + i + ".dat";
			var pathToLogFile = pathToRootDir + "/" + i + "/" + j + "/" + j
					+ ".log";
			var x_output = pathToRootDir + "/" + i + "/" + j + "/" + j
			+ ".xoutput";
			wfOut.processes.push(task("execute_case_" + i + "_" + j, functionName,
					pathToEllipsoidProgram, [ pathToDatFile, pathToLogFile,
							x_output ], [ "path_" + i + "_" + j ], [ "done_"
							+ i + "_" + j ])); // this executes ellipsoidy
		}
	}
	/*
	for (i = 0; i < iterationNumber; i++) { //name, functionName, ins, outs) {
		var doneTab =[]
		for (j = 0; j < numberOfBaseExecutions; j++) {
			doneTab.push("done_"+ i + "_" + j);
			doneTab.push("path_" + i + "_" + j);
		} // fulfils directories with dat files

		wfOut.processes.push(task("average_result_" + i, "createAverage", 
				pathToAveraging,["0"] ,doneTab, ["average_done_" + i ])); 
	}*/
	for (i = 0; i < iterationNumber; i++) {
		var xOutputTab = [];
		var doneTab =[]
		for (j = 0; j < numberOfBaseExecutions; j++) {
			doneTab.push("done_"+ i + "_" + j);
			doneTab.push("path_" + i + "_" + j);
			var x_output = pathToRootDir + "/" + i + "/" + j + "/" + j
			+ ".xoutput";
			xOutputTab.push(x_output)
		}
		xOutputTab.push(pathToRootDir + "/"+i+"/average.txt")
		wfOut.processes.push(task("average_result_" + i, functionName, 
				pathToAveraging,xOutputTab ,doneTab, ["average_done_" + i ])); 
	}
	for (i = 0; i < iterationNumber; i++) { //name, functionName, ins, outs) {

		wfOut.processes.push(functionExecution("pass_average_" + i, "passAverage", 
				["average_done_" + i,"path_" + i + "_" + "0" ], ["signal_permutation_average_" + i ])); 
	}
	
	
	
	var signalPermutationAverageTab =[]
	for (i = 0; i < iterationNumber; i++) { 
		signalPermutationAverageTab.push("signal_permutation_average_" + i);
		 
	}
	for (i = 0; i < iterationNumber; i++) { 
		signalPermutationAverageTab.push("signal_" + i +"_"+"0");
		 
	}
	signalPermutationAverageTab.push("rootPath");
	
	wfOut.processes.push(functionExecution("extraFunction" + i,
			"extraFunction", signalPermutationAverageTab, [
					"endMessage" ]));

	console.log(JSON.stringify(wfOut, null, 2)); // this gives on stdout
													// workflow
}

// here everything begins

// parser
if (!argv._[0] || !argv._[1] || !argv._[2] || !argv._[3]) {
	console.log("Usage: node generator.js path_to_root_dir path_to_ellipsoid_program number_of_same_data_executions path_to_averaging_program [command]");
	process.exit();
}

var command;

if (!argv._[4]) {
  command = "command";  
} else {
  command = argv._[4];
}
console.log(command);

// main function
createWf(command, argv._[0], argv._[1], argv._[2],argv._[3]);
