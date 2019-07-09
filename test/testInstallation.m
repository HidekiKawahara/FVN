%% Check for validity of installation

mypath = mfilename('fullpath');
myDirectory = mypath(1:strfind(mypath, 'test/testInstallation') - 1);

fs = 44100;
sigmaT = 0.01;
output = generateSparseFVNset(fs, sigmaT);
tmpCmd = ['load ' myDirectory 'src/fvnMin10ms.mat'];
eval(tmpCmd);

diffChk = output.fvnSet - fvnMin10ms;
if max(abs(diffChk(:))) == 0
    disp('Passed (generateStandardFVN, generateFVN6real, generateSparseFVNset)');
else
    disp('Saved data differs from the generated data.');
end
