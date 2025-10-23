@echo off
echo Activating Anaconda environment... 
call C:\ProgramData\anaconda3\Scripts\activate.bat C:\ProgramData\anaconda3 
echo Environment activated. Now running conda env update... 
conda env update -f \\sr-svw-914.unifr.ch\REPOSITORY\Pohl\sgg00425-master\_pystac_install-env.yml
echo Installation command completed. 
pause