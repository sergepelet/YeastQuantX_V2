Connection Between Matlab and Filemaker via JDBC:
Verify that the Database toolbox from Matlab is installed.
Copy fmjdbc.jar to /System/Library/Java/Extensions
or any folder in the javaclasspath of matlab



Filemaker scripts (Work only on Macs):

1.- Enable terminal to start matlab 
Edit the FMscript.command to add the full path of the matlab executable
 for instance: /Applications/MATLAB/MATLAB_R2009b.app/bin/matlab

or 

In the terminal copy the following line where 'xyz' is the matlab path (for instance: /Applications/MATLAB/MATLAB_R2009b.app)
echo 'export PATH=xyz/bin/matlab:$PATH' >> ~/.bash_profile


2.- YQlink.command

A) Edit the first line of the YQlink file to the actual path of the YeastQuantProgram/MatlabCode folder on your computer
B) copy YQlink.command in the ~/Library/FMconnect folder
