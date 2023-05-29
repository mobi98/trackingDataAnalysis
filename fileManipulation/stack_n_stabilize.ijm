// Choosing the directory containing the images to be stacked 
directory = getDirectory("Choose working directory");

folders = getFileList(directory);

for (i = 0; i < folders.length; i++){
	subdir = directory + folders[i]+'Ph/';
	folder_name = replace(folders[i], '/', "");
	output_dir = directory+folders[i]+folder_name+'_tl';
	stackNstabilize(subdir, output_dir);
	close(); 
}


function stackNstabilize(sd, od){
	
	run("Image Sequence...", "open=sd sort");
	run("Image Stabilizer", "transformation=Translation maximum_pyramid_levels=1 template_update_coefficient=0.90 maximum_iterations=200 error_tolerance=0.0000001");
	saveAs("Tiff", od);
}
