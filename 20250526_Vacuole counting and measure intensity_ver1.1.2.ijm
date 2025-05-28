/*20250506_Vacuole counting and measure intensity ver 1.1.1
 * Edited by Shao-Chun, Peggy, Hsu by 2025/3/12.
 * peggyschsu@ntu.edu.tw
 * - Draw cell contour and extract LC3 holes. 
 * - Get candida location from both DAPI and Candida channels.
 * - Get candida containing vacuoles (CCV).
 * - Measure LC3 intensity on CCV edges.

 * - ver 1.1.0 with following changes
 * ----1. limit the process file format with .tif, .lsm, .czi
 * ----2. Store all intermediate files in another folder created under the parent folder.
 * ----3. Report vacuole no. directly if no LC3 holes inside cells are detected.
 * ----4. Elevate the minimum size requirement for fungal DNA.

 * -ver 1.1.1 with following changes (editted by 2025/5/6)
 * ---1.Change the channel order from DAPI, Candida, actin,LC3 to DAPI, actin, Candida and LC3. 
 * ---2. Use actin channel to identify cell contour with area filter.
 * ---3. Candida was defined by the fungal DAPI surrounded by Candida GFP. Because some fungal DNA were not surrounded by Candida GFP,
 * ---these regions were defined by the dapi staining channel.
 * ---4. Cell no was counted by the actin channel to avoid the double nuclei over counting.
 
 * -ver1.1.2 (editted by 2025/5/27)
 * ---1. Count Candida cell wall inside of cells.
 * ---2. selectnone before filtering out Vacuole with low DAPI intensity 
 
 
 */
 
//Prepare for processing
run("Options...", "iterations=1 count=1 black");
//Prepare result table
    FilenameArray = newArray(0);
    CellnoArray = newArray(0);
    FungalDNAArray = newArray(0);
    FungalCellWallArray = newArray(0);
    VacuoleNoArray = newArray(0);
    LC3TotalInt = newArray(0);
    LC3IntAvg = newArray(0);
    Table.create("Vacoule analysis");
    Table.setColumn("File name", FilenameArray);
    Table.setColumn("Cell no.", CellnoArray);//actin channel
    Table.setColumn("Fungal DNA no.", FungalDNAArray);//dapi channel with size filter 0.3-10
    Table.setColumn("Fungal Cell Wall no.", FungalCellWallArray);//circles in Candida GFP channel    
    Table.setColumn("Vacuole no.", VacuoleNoArray);//OR(fungalDNA, fungalCellWall) 
    Table.setColumn("LC3 around vacuole", LC3TotalInt);//LC3 intensity sum with 4 pixels extension from each vacuole
    Table.setColumn("LC3 around vacuole/ vacuole", LC3IntAvg);//LC3 intensity average for each vacuole
//Loop for images
run("Set Measurements...", "area shape integrated redirect=None decimal=2");
dir1 = getDirectory("Folder of images");
list1 = getFileList(dir1);
dirP= File.getParent(dir1);
File.makeDirectory(dirP  + File.separator + "Intermediate and Result summary files");
path = dirP  + File.separator + "Intermediate and Result summary files";
for (j = 0; j < list1.length; j++) {
       	path1 = dir1 + list1[j];   
       	if(endsWith(list1[j], ".tif") || endsWith(list1[j], ".lsm")|| endsWith(list1[j], ".czi")){
        open(path1);
       	}
       	
//Prepare the image
    tifName = getTitle();
    if (endsWith(tifName, ".tif")) {
    	Name = replace(tifName, ".tif", "");
    }
     if (endsWith(tifName, ".lsm")) {
    	Name = replace(tifName, ".lsm", "");
    }
     if (endsWith(tifName, ".czi")) {
    	Name = replace(tifName, ".czi", "");
    }
    getDimensions(width, height, channels, slices, frames);  
    run("Set Measurements...", "area integrated redirect=None decimal=2");
    rename("Raw");
    run("Split Channels");
    selectImage("C1-Raw");
    rename("DAPI");
    selectImage("C3-Raw");
    rename("Candida");
    selectImage("C4-Raw");
    rename("LC3");
    selectImage("C2-Raw");
    rename("Actin");
/*    
//Draw cell outline to capture holes nearby cell membrane by LC3
  	selectImage("LC3");
  	run("Duplicate...", " ");
  	run("Gaussian Blur...", "sigma=3");
    setAutoThreshold("Mean dark no-reset");
    run("Analyze Particles...", "add");
    c = roiManager("count");
    a = newArray(0);
    for (i = 0; i < c ; i++) {
    	roiManager("select", i);
    	run("Convex Hull");  
    	roiManager("Add");
    	a = Array.concat(a,i);
    }
    roiManager("select", a);
    roiManager("delete");
    roiManager("Save", path + File.separator + Name+ "_Cell.zip");
  //Create cell mask
    newImage("CellMask", "8-bit black", width, height, 1);
    roiManager("Show All");
    setForegroundColor(255, 255, 255);
    roiManager("Fill");
    run("Erode");
    run("Erode");
    run("Divide...", "value=255");
    roiManager("reset");
*/
//Draw cell outline by actin
  	imageCalculator("Add create", "Actin","LC3");
	selectWindow("Result of Actin");
  	run("Gaussian Blur...", "sigma=4");
    setAutoThreshold("Mean dark no-reset");
    run("Analyze Particles...", "size=30-Infinity add");
    cellno = roiManager("count");
    roiManager("Save", path + File.separator + Name+ "_Cell.zip");
  //Create cell mask
    newImage("CellMask", "8-bit black", width, height, 1);
    setForegroundColor(255, 255, 255);
    for (i = 0; i < cellno; i++) {
    	selectWindow("Result of Actin");
    	roiManager("select", i);
    	run("Convex Hull");
    	selectWindow("CellMask");
    	run("Restore Selection");
    	run("Fill", "slice");
	}
	run("Select None");
    run("Divide...", "value=255");
    if (cellno>0) {
       	roiManager("reset");
       }

 
//Count cell no by DAPI stain
   	 selectImage("DAPI");
   	 run("Duplicate...", "title=Nucleus");
   	 selectImage("DAPI");
   	 run("Duplicate...", "title=DAPIRaw");
   	 selectWindow("Nucleus");
   	 run("Gaussian Blur...", "sigma=3");
     //imageCalculator("Multiply create", "Nucleus","CellMask");
     setAutoThreshold("Otsu dark");
     run("Analyze Particles...", "size=5-100 circularity=0.20-1.00 add");
     roiManager("Save", path + File.separator + Name +"_nucleus.zip");
     c = roiManager("count");
     /*selectWindow("Result of Nucleus");
     close();     */
     if (c>0) {
       	roiManager("reset");
       }

//Define location of Candida from smaller DAPI stain areas or labeled Candida
   //Generate Fungal DNA mask   
     selectWindow("DAPI");
     run("Gaussian Blur...", "sigma=2");
     setMinAndMax(0, 35);
     run("Apply LUT");
     run("Auto Local Threshold", "method=Phansalkar radius=2 parameter_1=0 parameter_2=0 white");
     run("Watershed");
     //run("Threshold...");
     run("Analyze Particles...", "size=0.30-10 circularity=[] show=Masks");
     imageCalculator("Multiply", "Mask of DAPI","CellMask");
     rename("FungalDNA");
     run("Grays");
     setAutoThreshold("Otsu dark");
     run("Analyze Particles...", "size=0.30-10.00 show=Masks add");
     FungalDNAno = roiManager("count");
     roiManager("Save", path + File.separator + Name + "_FungalDNA.zip");
     /*selectImage("FungalDNA");
     saveAs("tiff", path + File.separator + Name + "_FungalDNAMask.tif");
     selectWindow(Name + "_FungalDNAMask.tif");
     rename("FungalDNA");*/
     if (FungalDNAno>0) {
       	roiManager("reset");
     }
     resetThreshold;
 
  //Generate Candida region by Candida channel
	selectWindow("Candida");
	run("Median...", "radius=2");
	run("Find Maxima...", "prominence=8 light output=[Segmented Particles]");
	imageCalculator("Multiply", "Candida Segmented","CellMask");
	selectWindow("Candida Segmented");
	run("Analyze Particles...", "size=0-50 add");
	FungalGFPno = roiManager("count");
     if (FungalGFPno>0) {
       	roiManager("Save", path + File.separator + Name + "_FungalGFP.zip");
     }
    
    //Clear
  		selectWindow("Candida Segmented");
  		close();

  //Generate Candida location by filtering out FungalGFP without fungal DNA and keeping bigger region from either fungal DNA or fungal GFP 
  for (i = 0; i < FungalGFPno; i++) {
  		selectWindow("FungalDNA");
  		roiManager("Show None");
  		roiManager("select", i);
  		run("Measure");
  		selectWindow("Results");
  		fungalDNAInt = Table.get("RawIntDen", 0);
  		if (fungalDNAInt>500) {
  			roiManager("fill");
   		}
   		run("Clear Results");
  	}
	if (FungalGFPno>0) {
       	roiManager("reset");
    }
  //Filter Candida location outside of cells	
  	imageCalculator("Multiply", "FungalDNA","CellMask");
  //Filter Candida location with low DAPI signal	
    selectWindow("FungalDNA");
    run("Select None");
    run("Analyze Particles...", "size=0-50 add");
    c = roiManager("count");
    //print(c);
     for (i = 0; i < c; i++) {
       	selectWindow("DAPIRaw");
     	roiManager("select", i);
  		run("Measure");
  		selectWindow("Results");
  		DAPIInt = Table.get("RawIntDen", 0);
  		//print(i);
  		//print(DAPIInt);
  		if (DAPIInt<200) {
  			roiManager("select", i);
  			roiManager("delete");
  			c = c-1;
  			i = i-1;
  		}
   		run("Clear Results");
     }
     CandidaNo = roiManager("count");
     if (CandidaNo>0) {
       	roiManager("Save", path + File.separator + Name + "_Vacuole.zip");
     }
     /*newImage("CandidaLocation", "8-bit black", width, height, 1);
     roiManager("Show All");
     setForegroundColor(255, 255, 255);
     roiManager("fill");     
     saveAs("tiff", path + File.separator + Name + "_CandidaLocation.tif");*/

//Clear
	selectWindow("DAPI");
	close();
	selectImage("Actin");
	close();
	selectImage("Candida");
	close();
	selectImage("Result of Actin");
	close();
	selectImage("CellMask");
	close();
	selectImage("Nucleus");
	close();
	selectImage("DAPIRaw");
	close();
	selectImage("FungalDNA");
	close();
	selectImage("Mask of FungalDNA");
	close();
	
//Measure LC3 intensity around vacuoles
   CandidaNo =roiManager("count");
   if (CandidaNo>0) {
      newImage("Vacuole LC3", "8-bit black", width, height, 1);
      selectWindow("Vacuole LC3");
      roiManager("Set Color", "yellow");
	  roiManager("Set Line Width", 6);
      roiManager("show all without labels");
      run("Flatten");
      selectWindow("Vacuole LC3-1");
      run("8-bit");
      setAutoThreshold("Otsu dark");
      run("Convert to Mask");
      run("Dilate");
      run("Dilate");
      run("Divide...", "value=255");
      selectImage("Vacuole LC3");
      close();
      imageCalculator("Multiply create", "LC3","Vacuole LC3-1");
      selectWindow("Result of LC3");
      roiManager("Deselect");
      getRawStatistics(nPixels, mean, min, max, std, histogram);
      VacLC3Intensity = nPixels * mean;
      VacLC3IntensityAvg = nPixels * mean/CandidaNo;
      selectWindow("Result of LC3");
      saveAs("tiff", path + File.separator + Name + "_Extracted LC3 around vacuoles.tif");
    }
    else{
      VacLC3Intensity = 0;
      VacLC3IntensityAvg = 0;
    }
    //print(tifName + "      has   " + c + "   Vacuoles in   " + cellno + "   cells with average LC3 intensity " + VacLC3IntensityAvg + " /vacuole.");
    FilenameArray = Array.concat(FilenameArray,tifName);
    CellnoArray = Array.concat(CellnoArray,cellno);
    FungalDNAArray = Array.concat(FungalDNAArray, FungalDNAno);
    FungalCellWallArray = Array.concat(FungalCellWallArray, FungalGFPno);
    VacuoleNoArray = Array.concat(VacuoleNoArray,CandidaNo);
    LC3TotalInt = Array.concat(LC3TotalInt,VacLC3Intensity);
    LC3IntAvg = Array.concat(LC3IntAvg,VacLC3IntensityAvg);

//Clear
  	run("Close All");
  	roiManager("reset");
  	run("Clear Results");
  
//Generate result table
    selectWindow("Vacoule analysis");
    Table.setColumn("File name", FilenameArray);
    Table.setColumn("Cell no.", CellnoArray);
    Table.setColumn("Fungal DNA no.", FungalDNAArray);
    Table.setColumn("Fungal Cell Wall no.", FungalCellWallArray);
    Table.setColumn("Vacuole no.", VacuoleNoArray);
    Table.setColumn("LC3 around vacuole", LC3TotalInt);
    Table.setColumn("LC3 around vacuole/ vacuole", LC3IntAvg);
}
    saveAs("Results",  path + File.separator + "Vacuole analysis.csv" );