function [] = make_solnp()


% change your mkl path here
mkl_path = "C:\Program Files (x86)\Intel\oneAPI\mkl\2022.1.0";

% change your osqp path here
osqp_path = "E:\APP\osqp";

mkl_lib_path = mkl_path + "\lib\intel64";
osqp_lib_path = osqp_path + "\lib";


if exist(mkl_path) && exist(osqp_path)
    lib_path = [mkl_lib_path, osqp_lib_path];
    lib_path = join('-L' + lib_path);
    platform = convertCharsToStrings(computer('arch'));
else
    error('Please add your pathes to mkl and osqp\n');
end 


link = '';

if platform == "win64"
    fprintf('Linking MKL in Windows \n');
    link = [link, ' -lmkl_intel_lp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
elseif platform == "maci64"
    fprintf('Linking MKL in MacOS \n');
    link = [link, ' -lmkl_intel_ilp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
elseif platform == "glnxa64"
    fprintf('Linking MKL in Linux \n');
    link = [link, ' -lmkl_intel_ilp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
else
    error('Unsupported platform.\n');
end

link = [link, '-losqp.lib'];

mexfname = 'SOLNP';
psrc = '.\source';  

src_files = dir( [psrc '\*.c'] );
srclist = [];

for i = 1:length(src_files)
    srclist = [srclist,convertCharsToStrings(src_files(i).name)];
end

src = fullfile(psrc, srclist);

pinc = "include";


mkl_include = mkl_path + "\include";
osqp_include = osqp_path + "\include";

inc = [pinc,mkl_include,osqp_include];
inc = join('-I' + inc);

cmd = 'mex -O -output ' + join([mexfname, lib_path, src, inc, link]);

fprintf('%s\n',cmd);
eval(replace(cmd, 'Program Files (x86)', "'Program Files (x86)'"));

end