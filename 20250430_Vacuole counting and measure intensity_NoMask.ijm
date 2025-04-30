/*20250430_Vacuole counting and measure intensity
 * - Draw cell contour and extract LC3 holes. 
 * - Get candida location from both DAPI and Candida channels.
 * - Get candida containing vacuoles (CCV).
 * - Measure LC3 intensity on CCV edges.
 * - Compiled by Shao-Chun, Peggy, Hsu by 2025/3/12.
 * - peggyschsu@ntu.edu.tw
 * 
 * This edition skips the generation of intermediate image masks and ROIs for a more streamlined output.


 * 
 */
//Prepare result table
    FilenameArray = newArray(0);
    CellnoArray = newArray(0);
    VacuoleNoArray = newArray(0);
    LC3TotalInt = newArray(0);
    LC3IntAvg = newArray(0);
    Table.create("Vacoule analysis");
    Table.setColumn("File name", FilenameArray);
    Table.setColumn("Cell no.", CellnoArray);
    Table.setColumn("Vacuole no.", VacuoleNoArray);
    Table.setColumn("LC3 around vacuole", LC3TotalInt);
    Table.setColumn("LC3 around vacuole/ vacuole", LC3IntAvg);
//Loop for images
run("Set Measurements...", "area shape integrated redirect=None decimal=2");
dir1 = getDirectory("Folder of images");
list1 = getFileList(dir1);
for (j = 0; j < list1.length; j++) {
       	path1 = dir1 + list1[j];       
        open(path1);

//Prepare the image
    tifName = getTitle();
    path = getDirectory("image");
    getDimensions(width, height, channels, slices, frames);  
    run("Set Measurements...", "area integrated redirect=None decimal=2");
    rename("Raw");
    run("Split Channels");
    selectImage("C1-Raw");
    rename("DAPI");
    selectImage("C2-Raw");
    rename("Candida");
    selectImage("C4-Raw");
    rename("LC3");
    selectImage("C3-Raw");
    close();
//Draw cell outline to capture holes nearby cell membrane
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
    roiManager("Save", path + File.separator + "Cell.zip");
  //Create cell mask
    newImage("CellMask", "8-bit black", width, height, 1);
    roiManager("Show All");
    setForegroundColor(255, 255, 255);
    roiManager("Fill");
    run("Erode");
    run("Erode");
    run("Divide...", "value=255");
    roiManager("reset");
//Count cell no
   	 selectImage("DAPI");
   	 run("Duplicate...", "title=Nucleus");
     run("Gaussian Blur...", "sigma=3");
     imageCalculator("Multiply create", "Nucleus","CellMask");
     setAutoThreshold("Otsu dark");
     run("Analyze Particles...", "size=5-100 circularity=0.20-1.00 add");
     roiManager("Save", path + File.separator + "nucleus.zip");
     cellno = roiManager("count");
     //Clear
     selectWindow("Nucleus");
     close();
     selectWindow("Result of Nucleus");
     close();     
     if (cellno>0) {
       	roiManager("reset");
       }
//Get holes from LC3 channel
  //Get LC3 Hole
    selectWindow("LC3");
    run("Duplicate...", "title=LC3VacuoleIntensity");
    selectWindow("LC3");
    run("Enhance Contrast", "saturated=0.35");
    roiManager("Open", path + File.separator + "Cell.zip");
    roiManager("Show All without labels");
    run("Flatten");
    run("8-bit");
    c = roiManager("count");
       if (c>0) {
       	roiManager("reset");
       }
    // Init GPU
      run("CLIJ2 Macro Extensions", "cl_device=");
      Ext.CLIJ2_clear();
      input = "LC3-2";
      Ext.CLIJ2_push(input);
    // Tubeness
       sigma = 1.5;
       calibrationX = 2.0;
       calibrationY = 2.0;
       calibrationZ = 2.0;
       Ext.CLIJx_imageJ2Tubeness(input, output, sigma, calibrationX, calibrationY, calibrationZ);
       Ext.CLIJ2_release(input);
       Ext.CLIJ2_pull(output);
       Ext.CLIJ2_release(output);
       Ext.CLIJ2_clear();
     //Get LC3 Hole
       selectWindow(output);
       rename("Lateralized LC3");
       run("Find Maxima...", "prominence=1 light output=[Maxima Within Tolerance]");
       imageCalculator("Multiply create", "Lateralized LC3 Maxima","CellMask");
       selectWindow("Result of Lateralized LC3 Maxima");
       run("Open");
       setAutoThreshold("Mean dark no-reset");
       roiManager("Open", path + File.separator + "nucleus.zip");
       selectWindow("Result of Lateralized LC3 Maxima");
       roiManager("Show All without labels"); 
       setForegroundColor(0, 0, 0);
       roiManager("Fill");
       if (c>0) {
       	roiManager("reset");
       }
       
       run("Analyze Particles...", "size=80-1300 circularity=0.20-1.00 add");
       roiManager("Save", path + File.separator + "LC3Hole.zip");
       c = roiManager("count");
       if (c>0) {
       	roiManager("reset");
       }
    //Clear
      selectWindow("LC3-1");
      close();   
      selectImage("Lateralized LC3 Maxima");
      close();
      selectImage("Lateralized LC3");
      close();
      selectImage("LC3-2");
      close();
      selectImage("Result of Lateralized LC3 Maxima");
      close();
//Define location of Candida from smaller DAPI stain areas or labeled Candida
   //Generate Bac DNA mask   
     selectWindow("DAPI");
     run("Gaussian Blur...", "sigma=2");
     setMinAndMax(0, 42);
     run("Apply LUT");
     run("Auto Local Threshold", "method=Phansalkar radius=2 parameter_1=0 parameter_2=0 white");
     run("Despeckle");
     run("Watershed");
     //run("Threshold...");
     run("Analyze Particles...", "size=0.050-10.00 add");
     roiManager("Save", path + File.separator + "bacDNA.zip");
     newImage("BacDNA", "8-bit black", width, height, 1);
     roiManager("Show All with labels");
     setForegroundColor(255, 255, 255);
     roiManager("Fill");
     c = roiManager("count");
     if (c>0) {
       	roiManager("reset");
     }
   //Generate Candida location from DAPI and Candida channel
     selectImage("Candida");
     run("Gaussian Blur...", "sigma=1");
     setAutoThreshold("MaxEntropy dark no-reset");
     run("Convert to Mask");
     imageCalculator("Add create", "Candida","BacDNA");
     selectImage("Result of Candida");
     rename("CandidaLocation");
   //Remove the cell nucleus from CandidaLocation mask
     selectWindow("CandidaLocation");
     roiManager("Open", path + File.separator + "nucleus.zip");
     selectWindow("CandidaLocation");
     roiManager("Show All without labels"); 
     setForegroundColor(0, 0, 0);
     roiManager("Fill");
     saveAs("tiff", path + File.separator + "CandidaLocation.tif");
   //Document 
     selectImage("Candida");
     saveAs("tiff", path + File.separator + "CandidaMask.tif");
     close();
     selectImage("BacDNA");
     saveAs("tiff", path + File.separator + "BacIMask.tif");
     close();
     selectImage("DAPI");
     close();
     c = roiManager("count");
       if (c>0) {
       	roiManager("reset");
       }
//Justify LC3 surrounded Candida
  roiManager("Open", path + File.separator + "LC3Hole.zip");
  roiManager("Measure");
  c =roiManager("count");
  nonVC = newArray(0);
  for (i = 0; i < c; i++) {
  	  selectWindow("Results");
  	  tester = getResult("RawIntDen", i);
  	  if (tester < 511) {
  		  nonVC = Array.concat(nonVC,i);
  	  }
    }
    //Array.print(nonVC);
    if (nonVC.length >0) {
    	 roiManager("select", nonVC);
         roiManager("delete");
         c =roiManager("count");
         if (c>0) {
         	roiManager("Save", path + File.separator + "Vacuole.zip");
         	newImage("Vacuole LC3", "8-bit black", width, height, 1);
         	selectWindow("Vacuole LC3");
         	roiManager("show all without labels");
         	roiManager("Set Line Width", 6);
         	run("Flatten");
         	selectWindow("Vacuole LC3-1");
         	run("8-bit");
         	run("Divide...", "value=170");
         	selectImage("Vacuole LC3");
            close();
            imageCalculator("Multiply create", "LC3VacuoleIntensity","Vacuole LC3-1");
            selectWindow("Result of LC3VacuoleIntensity");
            roiManager("Deselect");
            getRawStatistics(nPixels, mean, min, max, std, histogram);
            VacLC3Intensity = nPixels * mean;
            VacLC3IntensityAvg = nPixels * mean/c;
            selectWindow("Result of LC3VacuoleIntensity");
            saveAs("tiff", path + File.separator + "Extracted LC3 around vacuoles.tif");
         }
         else{
         	VacLC3Intensity = 0;
         	VacLC3IntensityAvg = 0;
         }
    }
    print(tifName + "      has   " + c + "   Vacuoles in   " + cellno + "   cells with average LC3 intensity " + VacLC3IntensityAvg + " /vacuole.");
    FilenameArray = Array.concat(FilenameArray,tifName);
    CellnoArray = Array.concat(CellnoArray,cellno);
    VacuoleNoArray = Array.concat(VacuoleNoArray,c);
    LC3TotalInt = Array.concat(LC3TotalInt,VacLC3Intensity);
    LC3IntAvg = Array.concat(LC3IntAvg,VacLC3IntensityAvg);
//Clear
  	run("Close All");
  	roiManager("reset");
  	run("Clear Results");
  	
}
//Generate result table
    selectWindow("Vacoule analysis");
    Table.setColumn("File name", FilenameArray);
    Table.setColumn("Cell no.", CellnoArray);
    Table.setColumn("Vacuole no.", VacuoleNoArray);
    Table.setColumn("LC3 around vacuole", LC3TotalInt);
    Table.setColumn("LC3 around vacuole/ vacuole", LC3IntAvg);
    saveAs("Results",  path + File.separator + "Vacuole analysis.csv" );