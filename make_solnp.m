function [] = make_solnp()

mexfname = 'SOLNP_plus';

% change your mkl path here
% mkl_path = "F:\MKL\mkl\2021.3.0";
mkl_path = '/nfsshare/software/deng/install/intel/mkl';

% change your osqp path here
% osqp_path = "F:\GitHub\SOLNP\osqp";
osqp_path = '/nfsshare/home/liujinsong/github/osqp';

% mkl_lib_path = mkl_path + "\lib\intel64";
% osqp_lisb_path = osqp_path + "\lib";

mkl_lib_path = fullfile(mkl_path, 'lib', 'intel64');
osqp_lib_path = fullfile(osqp_path, 'lib');



if exist(mkl_path) && exist(osqp_path)
    lib_path = [convertCharsToStrings(mkl_lib_path), convertCharsToStrings(osqp_lib_path)];
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
    link = [link, ' -lmkl_intel_lp64',...
        ' -lmkl_core', ' -lmkl_sequential '];
elseif platform == "glnxa64"
    fprintf('Linking MKL in Linux \n');
    link = [link, ' -lmkl_intel_lp64',...
        ' -lmkl_core', ' -lmkl_sequential '];

    sys_cmd = ['gcc -fPIC -I', matlabroot, '/extern/include ' , '-I./include -I', mkl_path, '/include -I', osqp_path, '/include ',...
                '-shared -o ', mexfname, '.mexa64 ', './source/*.c ', '-Wl,--start-group ', mkl_lib_path, '/libmkl_core.a ',...
                mkl_lib_path, '/libmkl_intel_lp64.a ', mkl_lib_path, '/libmkl_sequential.a ', mkl_lib_path, '/libmkl_intel_thread.a ',...
                osqp_lib_path, '/libosqp.a -Wl,--end-group -L', matlabroot, '/bin/glnxa64 -ldl -lpthread -lmx -lmex -lmat -lm'];

    system(sys_cmd);
    return
else
    error('Unsupported platform.\n');
end

link = [link, '-losqp'];

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

cmd = 'mex -g -output ' + join([mexfname, lib_path, src, inc, link]);

fprintf('%s\n',cmd);
eval(replace(cmd, 'Program Files (x86)', "'Program Files (x86)'"));

end