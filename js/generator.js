var argv = require('optimist').argv;

function createBasicSignals(wfOut, pathToRootDir, pathToEllipsoidProgram){
	var signals = []
	signals.push("rootPath")
	signals.push("ellipsoidPath")
	signals.push("datDirs")
	
	for (i = 0; i < signals.length; i++) {
		wfOut.signals.push({
			name : signals[i]
		});
	}
	
	wfOut.signals[0].data = [ pathToRootDir ]
	wfOut.signals[1].data = [ pathToEllipsoidProgram ]
	wfOut.signals[2].data = [ ]
}
function createUsualSignals(howMuch, prefix_name, wfOut){
	
	for (i = 0; i < howMuch; i++) {
		wfOut.signals.push({
			name : prefix_name+i,
			data : []
		});
	}
}



function createSignal(wfOut, index, ifCoord,numberOfParticles,
		numberOfSpecies,forceScalingFactor,
		rotationScalingFactor,diameterIncreasingFactor,
		cellsX,cellsY,
		cellsZ,constractionRate,
		initialPackingDensity,maximalNumberOfSteps,
		isConstantVolume,numberOfParts,
		diameterOfParts1,diameterOfParts2,
		diameterOfParts3,numberOfLinesPerPage,
		numberOfStepsBetweenPrintouts,
		numberOfStepsBetweenCoord,
		numberOfStepsBetweenRotations){
	
	wfOut.signals.push({name : "signal"+index});
	index+=3
	wfOut.signals[index].ifCoord =  [ifCoord] 
	wfOut.signals[index].numberOfParticles =  [numberOfParticles] 
	wfOut.signals[index].numberOfSpecies =  [numberOfSpecies] 
	wfOut.signals[index].forceScalingFactor = [ forceScalingFactor] 
	wfOut.signals[index].rotationScalingFactor =  [rotationScalingFactor] 
	wfOut.signals[index].diameterIncreasingFactor =  [diameterIncreasingFactor] 
	wfOut.signals[index].cellsX =  [cellsX] 
	wfOut.signals[index].cellsY =  [cellsY] 
	wfOut.signals[index].cellsZ =  [cellsZ] 
	wfOut.signals[index].constractionRate = [ constractionRate] 
	wfOut.signals[index].initialPackingDensity =[  initialPackingDensity] 
	wfOut.signals[index].maximalNumberOfSteps = [ maximalNumberOfSteps] 
	wfOut.signals[index].isConstantVolume = [isConstantVolume] 
	wfOut.signals[index].numberOfParts =  [numberOfParts] 
	wfOut.signals[index].diameterOfParts1 =  [diameterOfParts1] 
	wfOut.signals[index].diameterOfParts2 =  [diameterOfParts2 ]
	wfOut.signals[index].diameterOfParts3 =  [diameterOfParts3] 
	wfOut.signals[index].numberOfStepsBetweenPrintouts =  [numberOfStepsBetweenPrintouts] 
	wfOut.signals[index].numberOfStepsBetweenCoord =  [numberOfStepsBetweenCoord ]
	wfOut.signals[index].numberOfStepsBetweenRotations =  [numberOfStepsBetweenRotations] 
	
}

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



function functionExecution(name, functionName, ins, outs) {
	return {
		"name" : name,
		"function" : functionName,
		"type" : "dataflow",
		"ins" : ins,
		"outs" : outs
	}
}
function createWf(functionName, pathToRootDir, pathToEllipsoidProgram) {
	//main tuple
	var wfOut = {
		processes : [],
		signals : [],
		ins : [],
		outs : [ 0 ]
	};
	
	// data variables
	var ifCoord = ["1"];
	var numberOfParticles = ["1000"];
	var numberOfSpecies = ["1"];
    var forceScalingFactor =  ["0.1"];
	var rotationScalingFactor = ["3.00"];
	var diameterIncreasingFactor = ["0.01"];
    var cellsX = ["20","40"];
    var cellsY = ["20"];
	var cellsZ = ["20"];
    var constractionRate = ["102400"];
	var initialPackingDensity = ["0.8"];
	var maximalNumberOfSteps = ["10000000"];
	var isConstantVolume = ["1"];
	var numberOfParts = ["1000"];
	var diameterOfParts1 = ["1.3"];
    var	diameterOfParts2 = ["1.0"];
	var diameterOfParts3 = ["1.0"];
	var numberOfLinesPerPage = ["56"];
	var numberOfStepsBetweenPrintouts = ["100"];
	var numberOfStepsBetweenCoord = ["1000000"];
    var	numberOfStepsBetweenRotations = ["1"];
	
    var lengths = [ifCoord.length,numberOfParticles.length,numberOfSpecies.length,
           		forceScalingFactor.length, rotationScalingFactor.length, 
           		diameterIncreasingFactor.length, cellsX.length, cellsY.length, cellsZ.length, 
           		constractionRate.length, initialPackingDensity.length,maximalNumberOfSteps.length,
           		isConstantVolume.length, numberOfParts.length, diameterOfParts1.length, 
           		diameterOfParts2.length, diameterOfParts3.length, numberOfLinesPerPage.length, 
           		numberOfStepsBetweenPrintouts.length,numberOfStepsBetweenCoord.length, 
           		numberOfStepsBetweenRotations.length];
    
    var counters = [lengths.length];
    var product = 1;
	for (i = 0; i < lengths.length; i ++) {
	    counters[i] = lengths[i];
	    product *= lengths[i];
	}
	var pathToRoot = String(pathToRootDir)
	var iterationNumber = 0 
	//ins
	wfOut.ins.push("rootPath");
	
	//signals
	
	createBasicSignals(wfOut,pathToRootDir, pathToEllipsoidProgram);
	//loop below creates signal with data for alll permutation
	for(i = 0; i < product; i++){
		var datDirName = pathToRoot+'/'+i;
		//console.log("This is path to dat dir name " + datDirName);
		var index = 0;
		createSignal(wfOut,iterationNumber,ifCoord[counters[index++] - 1],numberOfParticles[counters[index++] - 1],
			numberOfSpecies[counters[index++] - 1],forceScalingFactor[counters[index++] - 1],
			rotationScalingFactor[counters[index++] - 1],diameterIncreasingFactor[counters[index++] - 1],
			cellsX[counters[index++] - 1],cellsY[counters[index++] - 1],
			cellsZ[counters[index++] - 1],constractionRate[counters[index++] - 1],
			initialPackingDensity[counters[index++] - 1],maximalNumberOfSteps[counters[index++] - 1],
			isConstantVolume[counters[index++] - 1],numberOfParts[counters[index++] - 1],
			diameterOfParts1[counters[index++] - 1],diameterOfParts2[counters[index++] - 1],
			diameterOfParts3[counters[index++] - 1],numberOfLinesPerPage[counters[index++] - 1],
			numberOfStepsBetweenPrintouts[counters[index++] - 1],
			numberOfStepsBetweenCoord[counters[index++] - 1],
			numberOfStepsBetweenRotations[counters[index++] - 1]
		);
		iterationNumber+=1
		
		counters[lengths.length - 1]--;
		for(j = lengths.length -1; j > 0; j--){
			if(counters[j] == 0){
				counters[j] = lengths[j];
				counters[j-1]--;
			}
		}
	}
	createUsualSignals(iterationNumber,"path",wfOut);
	createUsualSignals(iterationNumber,"done",wfOut);
	var signalNames = []
	for(i = 0; i < iterationNumber; i++){
		signalNames.push("signal"+i) 
	}

	//processes
	wfOut.processes.push(functionExecution("mkdir", "createRootDir", 
			[0,"rootPath"], signalNames));  //create dat
	//wfOut.processes.push(task("createWorkingDirectory", functionName, 
		//	"mkdir", pathToRootDir,[], [])) //mkdir
	for(i = 0; i < iterationNumber; i++){ //(name, functionName, ins, outs)
		wfOut.processes.push(functionExecution("create_dat"+i, "createDat", 
				["signal"+i,"rootPath"], ["path"+i,"datDirs"]));  //create dats
	}
	for(i = 0; i < iterationNumber; i++){ //(name, functionName, executable, args, ins, outs) 
		var pathToDatFile = pathToRootDir + "/" + i + "/" + i+".dat";
		var pathToLogFile = pathToRootDir + "/" + i + "/" + i+".log";
		var x_output = "1"
		wfOut.processes.push(task("execute_case_"+i,"command", pathToEllipsoidProgram, 
				[pathToDatFile,pathToLogFile,x_output],
				["path"+i], ["done"+i]));  //create dats
	}
	//wfOut.processes.push(task("amqpCommand", functionName, 
		//	"mkdir", pathToRootDir)); //execute queries
	
	//return workflow
	console.log(JSON.stringify(wfOut, null, 2));
}


if (!argv._[0] || !argv._[1]) {
	console.log("Usage: node generator.js path_to_root_dir path_to_ellipsoid_program");
	process.exit();
}

createWf("command", argv._[0], argv._[1]);

