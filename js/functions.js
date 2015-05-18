var argv = require('optimist').argv;
var mkdirp = require('mkdirp');
var fs = require('fs');

function createDatString(ifCoord, numberOfParticles, numberOfSpecies,
		forceScalingFactor, rotationScalingFactor, diameterIncreasingFactor,
		cellsX, cellsY, cellsZ, constractionRate, initialPackingDensity,
		maximalNumberOfSteps, isConstantVolume, numberOfParts,
		diameterOfParts1, diameterOfParts2, diameterOfParts3,
		numberOfLinesPerPage, numberOfStepsBetweenPrintouts,
		numberOfStepsBetweenCoord, numberOfStepsBetweenRotations) {
	var string = '1.1  BASIC CONTROL DATA\n'
			+ '=true if coordinates are to be saved    (         Lreslt_x)      '
			+ ifCoord
			+ '\n'
			+ '2.1  PRINCIPAL SYSTEM PARAMETERS\n'
			+ 'Number of particles                     (         No_parts)         '
			+ numberOfParticles
			+ '\n'
			+ 'Number of species                       (       No_species)         '
			+ numberOfSpecies
			+ '\n'
			+ 'Force scaling factor                    (          Epsilon)         '
			+ forceScalingFactor
			+ '\n'
			+ 'Rotation scaling factor                 (          Eps_rot)         '
			+ rotationScalingFactor
			+ '\n'
			+ 'Diameter increasing factor              (        Diam_incr)         '
			+ diameterIncreasingFactor
			+ '\n'
			+ 'Number of cells in x direction          (       No_cells_x)         '
			+ cellsX
			+ '\n'
			+ 'Number of cells in y direction          (       No_cells_y)         '
			+ cellsY
			+ '\n'
			+ 'Number of cells in z direction          (       No_cells_z)         '
			+ cellsZ
			+ '\n'
			+ 'Contraction rate of outer diameter      (             Ntau)         '
			+ constractionRate
			+ '\n'
			+ 'Initial packing density                 (            Pnom0)         '
			+ initialPackingDensity
			+ '\n'
			+ 'Maximal number of steps                 (        Max_steps)     	 '
			+ maximalNumberOfSteps
			+ '\n'
			+ '= true if constant volume run          (          Leq_vol)         '
			+ isConstantVolume
			+ '\n'
			+ 'Number of parts (species 0)             (           Number)         '
			+ numberOfParts
			+ '\n'
			+ 'Diameter of parts (species 0)           (         Diameter)         '
			+ diameterOfParts1
			+ '\n'
			+ 'Diameter of parts (species 0)           (         Diameter)         '
			+ diameterOfParts2
			+ '\n'
			+ 'Diameter of parts (species 0)           (         Diameter)         '
			+ diameterOfParts3
			+ '\n'
			+ '4.1  PRINTOUT CONTROL\n'
			+ 'Number of lines per page                (        Npage_len)         '
			+ numberOfLinesPerPage
			+ '\n'
			+ 'number of steps between printouts       (      Nprint_step)         '
			+ numberOfStepsBetweenPrintouts
			+ '\n'
			+ 'number of steps between coord. storage  (       Nrslt_step)         '
			+ numberOfStepsBetweenCoord
			+ '\n'
			+ 'number of steps between rotations       (        Nrot_step)         '
			+ numberOfStepsBetweenRotations;
	return string;
}
function createDats(ins, outs, config, cb) {
	var ifCoord = [ "1" ];
	var numberOfParticles = [ "1000" ];
	var numberOfSpecies = [ "1" ];
	var forceScalingFactor = [ "0.1" ];
	var rotationScalingFactor = [ "3.00" ];
	var diameterIncreasingFactor = [ "0.01" ];
	var cellsX = [ "20", "40" ];
	var cellsY = [ "20", "40" ];
	var cellsZ = [ "20", "40" ];
	var constractionRate = [ "102400" ];
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
	var pathToRoot = String(ins.rootPath.data[0])
	mkdirp(pathToRoot, function(err) {
	});
	outs.datDirs.data = [];
	for (i = 0; i < product; i++) {
		var datDirName = pathToRoot + '/' + i;
		mkdirp(datDirName, function(err) {
		});
		var index = 0;
		var myFile = createDat(ifCoord[counters[index++] - 1],
				numberOfParticles[counters[index++] - 1],
				numberOfSpecies[counters[index++] - 1],
				forceScalingFactor[counters[index++] - 1],
				rotationScalingFactor[counters[index++] - 1],
				diameterIncreasingFactor[counters[index++] - 1],
				cellsX[counters[index++] - 1], cellsY[counters[index++] - 1],
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
		var datFilePath = datDirName + '/' + i + '.dat';
		var checkFile = -1;
		fs.writeFile(datFilePath, myFile, function(err) {
			if (err) {
				console.log("Couldnt create " + datFilePath + " file");
			} else {

			}
		});
		outs.datDirs.data.push(datFilePath);

		counters[lengths.length - 1]--;
		for (j = lengths.length - 1; j > 0; j--) {
			if (counters[j] == 0) {
				counters[j] = lengths[j];
				counters[j - 1]--;
			}
		}

	}
	cb(null, outs);
}
function createRootDir(ins, outs, config, cb) {
	var pathToRoot = String(ins.rootPath.data);
	mkdirp(pathToRoot, function(err) {
	});
	cb(null, outs);
}

function passAverage(ins, outs, config, cb) {
	var average
	var tab = String(ins[0].data).split("/");
	var i = 0;
	var path_averaged = "";
	for (i = 1; i < tab.length - 1; i++) {
		path_averaged = path_averaged + "/" + tab[i];
	}
	
	var fileAveraged = path_averaged + "/average.txt";
	average =fs.readFileSync(fileAveraged,"utf8" )
	outs[0].data = [];
	outs[0].data.push(average);

	cb(null, outs);

}
function extraFunction(ins, outs, config, cb) {

	var sumFilePath = ins[0].data + "/summary.csv";

	var columnNames = "ifCoord,numberOfParticles,numberOfSpecies,forceScalingFactor,"
			+ "rotationScalingFactor,diameterIncreasingFactor,"
			+ "cellsX,cellsY,cellsZ,constractionRate,initialPackingDensity,"
			+ "maximalNumberOfSteps,isConstantVolume,numberOfParts,"
			+ "diameterOfParts1,diameterOfParts2,diameterOfParts3,"
			+ "numberOfLinesPerPage,numberOfStepsBetweenPrintouts,numberOfStepsBetweenCoord,"
			+ "numberOfStepsBetweenRotations,averaged" + "\n";
	fs.writeFile(sumFilePath, columnNames, function(err) {
		if (err) {
			console.log("Couldnt create " + columnNames + " file");
		} else {
		}
	});
	for (i = 1; i < ins.length / 2; i++) {
		var columnValues = ins[i].ifCoord + "," + ins[i].numberOfParticles
				+ "," + ins[i].numberOfSpecies + ","
				+ ins[i].forceScalingFactor + ","
				+ ins[i].rotationScalingFactor + ","
				+ ins[i].diameterIncreasingFactor + "," + ins[i].cellsX + ","
				+ ins[i].cellsY + "," + ins[i].cellsZ + ","
				+ ins[i].constractionRate + "," + ins[i].initialPackingDensity
				+ "," + ins[i].maximalNumberOfSteps + ","
				+ ins[i].isConstantVolume + "," + ins[i].numberOfParts + ","
				+ ins[i].diameterOfParts1 + "," + ins[i].diameterOfParts2 + ","
				+ ins[i].diameterOfParts3 + "," + ins[i].numberOfLinesPerPage
				+ "," + ins[i].numberOfStepsBetweenPrintouts + ","
				+ ins[i].numberOfStepsBetweenCoord + ","
				+ ins[i].numberOfStepsBetweenRotations + ","
				+ ins[(i + ((ins.length - 1) / 2))].data[0] + "\n";

		fs.appendFile(sumFilePath, columnValues, function(err) {
		});

	}

	outs[0].data = [];
	outs[0].data.push("This is end of hyperflw workflow execution");
	console.log("This is end of hyperflw workflow execution");
	// proces.exit();
	cb(null, outs);
}
function createDat(ins, outs, config, cb) {

	var index_iteration = (ins[1].name.split("signal"))[1].split("_");
	var pathToRoot = String(ins.rootPath.data[0])
	var datDirPath = pathToRoot + '/' + index_iteration[1];
	var iterationPath = datDirPath + '/' + index_iteration[2];

	if (!fs.existsSync(datDirPath)) {
		mkdirp(datDirPath, function(err) {
		});
	}
	mkdirp(iterationPath, function(err) {
	});

	var myFile = createDatString(ins[1].ifCoord, ins[1].numberOfParticles,
			ins[1].numberOfSpecies, ins[1].forceScalingFactor,
			ins[1].rotationScalingFactor, ins[1].diameterIncreasingFactor,
			ins[1].cellsX, ins[1].cellsY, ins[1].cellsZ,
			ins[1].constractionRate, ins[1].initialPackingDensity,
			ins[1].maximalNumberOfSteps, ins[1].isConstantVolume,
			ins[1].numberOfParts, ins[1].diameterOfParts1,
			ins[1].diameterOfParts2, ins[1].diameterOfParts3,
			ins[1].numberOfLinesPerPage, ins[1].numberOfStepsBetweenPrintouts,
			ins[1].numberOfStepsBetweenCoord,
			ins[1].numberOfStepsBetweenRotations);
	var datFilePath = datDirPath + '/' + index_iteration[1] + '.dat';
	fs.writeFile(datFilePath, myFile, function(err) {
		if (err) {
			console.log("Couldnt create " + datFilePath + " file");
		} else {
		}
	});

	outs[0].data = [];
	outs[0].data.push(iterationPath);
	outs[1].data = [];
	outs[1].data.push([ datFilePath,
			"signal_" + index_iteration[1] + "_" + index_iteration[2] ]);
	cb(null, outs);
}

exports.createDat = createDat;
exports.createRootDir = createRootDir;
exports.passAverage = passAverage;
exports.extraFunction = extraFunction;